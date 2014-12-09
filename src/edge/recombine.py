from django.db import transaction
from edge.blast import blast_genome
from edge.models import Fragment, Operation
from Bio.Seq import Seq


class RecombinationRegion(object):

    def __init__(self, fragment_id, fragment_name,
                 start, end, sequence, cassette_reversed, front_arm, back_arm):
        self.fragment_id = fragment_id
        self.fragment_name = fragment_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.cassette_reversed = cassette_reversed
        self.front_arm = front_arm
        self.back_arm = back_arm

    def to_dict(self):
        return self.__dict__


def compute_swap_region_from_results(front_arm_sequence, front_arm_blastres,
                                     back_arm_sequence, back_arm_blastres):
    """
    Computes a region on fragment flanked by two arms from blast results.
    """

    MIN_IDENTITIES = 0.9
    MAX_MISSING_BP = 3

    # arms must be on the same fragment
    if front_arm_blastres.fragment_id != back_arm_blastres.fragment_id:
        return None

    # arms must align to same strand
    front_arm_strand = front_arm_blastres.strand()
    back_arm_strand = back_arm_blastres.strand()
    if front_arm_strand != back_arm_strand:
        return None

    # arms must align nicely
    if front_arm_blastres.identity_ratio() < MIN_IDENTITIES or\
       back_arm_blastres.identity_ratio() < MIN_IDENTITIES or\
       front_arm_blastres.alignment_length() < len(front_arm_sequence)-MAX_MISSING_BP or\
       back_arm_blastres.alignment_length() < len(back_arm_sequence)-MAX_MISSING_BP:
        return None

    # must align in the right orientation
    if front_arm_blastres.subject_start > front_arm_blastres.subject_end:
        # aligns to antisense strand
        if back_arm_blastres.subject_start > front_arm_blastres.subject_end:
            return None
        region_start = back_arm_blastres.subject_end
        region_end = front_arm_blastres.subject_start
        is_reversed = True
    else:
        # aligns to sense strand
        if back_arm_blastres.subject_start < front_arm_blastres.subject_end:
            return None
        region_start = front_arm_blastres.subject_start
        region_end = back_arm_blastres.subject_end
        is_reversed = False

    # get sequence between arms, including arm
    fragment = front_arm_blastres.fragment.indexed_fragment()
    region = fragment.get_sequence(bp_lo=region_start, bp_hi=region_end)

    return RecombinationRegion(fragment.id,
                               fragment.name,
                               region_start,
                               region_end,
                               region,
                               is_reversed,
                               front_arm_sequence,
                               back_arm_sequence)


def find_possible_swap_regions(genome, front_arm_sequence, back_arm_sequence):
    """
    Computes all possible recombined region from arm sequences.
    """

    front_arm_results = blast_genome(genome, 'blastn', front_arm_sequence)
    back_arm_results = blast_genome(genome, 'blastn', back_arm_sequence)

    regions = []
    for a_res in front_arm_results:
        for b_res in back_arm_results:
            region = compute_swap_region_from_results(front_arm_sequence, a_res,
                                                      back_arm_sequence, b_res)
            if region is not None:
                regions.append(region)

    return regions


def find_swap_region(genome, cassette, min_homology_arm_length):
    """
    Find a region on genome that can be recombined out using the cassette.
    Returns homology arms and possible regions.
    """

    if len(cassette) < 2*min_homology_arm_length:
        return None

    front_arm = cassette[0:min_homology_arm_length]
    back_arm = cassette[-min_homology_arm_length:]

    return find_possible_swap_regions(genome, front_arm, back_arm)


@transaction.atomic()
def recombine(genome, cassette, min_homology_arm_length,
              genome_name=None, cassette_name=None, notes=None):
    cassette = str(Seq(cassette))  # clean the sequence

    regions = find_swap_region(genome, cassette, min_homology_arm_length)
    if regions is None or len(regions) != 1:
        return None

    new_genome = genome.update()

    new_fragment_id = None
    with new_genome.update_fragment_by_fragment_id(regions[0].fragment_id) as f:
        new_region = cassette
        if regions[0].cassette_reversed:
            new_region = str(Seq(new_region).reverse_complement())
        f.replace_bases(regions[0].start, regions[0].end-regions[0].start+1, new_region)
        new_fragment_id = f.id

    cassette_name = 'Recombination cassette' if cassette_name is None else cassette_name
    with new_genome.annotate_fragment_by_fragment_id(new_fragment_id) as f:
        f.annotate(regions[0].start, regions[0].start+len(new_region)-1,
                   cassette_name, 'feature', 1)

    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s recombined %d-%d with %d bps" % (genome.name,
                                                           regions[0].start,
                                                           regions[0].end,
                                                           len(cassette))
    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    fragment = Fragment.create_with_sequence(name=cassette_name,
                                             circular=False,
                                             sequence=cassette)
    op = Operation(type=Operation.RECOMBINATION[0], fragment=fragment)
    op.save()
    new_genome.operations.add(op)

    return new_genome
