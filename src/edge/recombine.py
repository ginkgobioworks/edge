import json
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

    fragment = front_arm_blastres.fragment.indexed_fragment()

    # must align in the right orientation
    if front_arm_blastres.strand() < 0:
        # aligns to antisense strand
        if back_arm_blastres.subject_start > front_arm_blastres.subject_end and\
           fragment.circular is False:
            return None
        region_start = back_arm_blastres.subject_end
        region_end = front_arm_blastres.subject_start
        is_reversed = True
    else:
        # aligns to sense strand
        if back_arm_blastres.subject_start < front_arm_blastres.subject_end and\
           fragment.circular is False:
            return None
        region_start = front_arm_blastres.subject_start
        region_end = back_arm_blastres.subject_end
        is_reversed = False

    # get sequence between arms, including arm
    region_start = fragment.circ_bp(region_start)
    region_end = fragment.circ_bp(region_end)

    if region_end < region_start:
        assert fragment.circular is True
        region_p1 = fragment.get_sequence(bp_lo=region_start)
        region_p2 = fragment.get_sequence(bp_hi=region_end)
        region = region_p1+region_p2

    else:
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


def recombine_region(genome, region, cassette, min_homology_arm_length, cassette_name,
                     op, need_new_fragment):
    """
    Recombines on a given region. Returns recombination cassette location, how
    many base pairs removed, how many base pairs added, new fragment id.
    """

    region_start = region.start
    region_end = region.end

    new_fragment_id = None
    with genome.update_fragment_by_fragment_id(region.fragment_id,
                                               new_fragment=need_new_fragment) as f:
        if region.cassette_reversed:
            cassette = str(Seq(cassette).reverse_complement())

        replaced = 0
        if region_start < region_end:
            f.replace_bases(region_start, region_end-region_start+1, cassette)
            replaced = region_end-region_start+1
        else:
            assert f.circular is True
            f.replace_bases(region_start, f.length-region_start+1, cassette)
            replaced = f.length-region_start+1
            f.remove_bases(1, region_end)
            replaced += region_end
            # adjust region_start after removing sequence at the start
            region_start -= region_end

        new_fragment_id = f.id

    with genome.annotate_fragment_by_fragment_id(new_fragment_id) as f:
        # region_start is already adjusted
        f.annotate(region_start, region_start+len(cassette)-1,
                   cassette_name, 'operation', 1, operation=op)

    return (region_start, replaced, len(cassette), new_fragment_id)


def shift_regions(regions, fragment_id, start, replaced, added, new_fragment_id):
    """
    For each region in regions list on the specified fragment, shift its
    coordinates based on number of base pairs replaced and added. Expects list
    of regions to be sorted by starting position and 'start' to be less than
    starting position of first region on the list.
    """

    if len(regions) == 0:
        return []

    assert start <= regions[0].start

    def updated_region(region):
        if region.start <= start+replaced-1:  # overlapping
            return None
        region.fragment_id = new_fragment_id
        if region.start < region.end:
            region.start = region.start-replaced+added
            region.end = region.end-replaced+added
        else:
            # last replacement cannot be across circular boundary, since we are
            # replacing in order of starting position
            region.start = region.start-replaced+added
        return region

    regions = [region if region.fragment_id != fragment_id else updated_region(region)
               for region in regions]
    return [r for r in regions if r is not None]


@transaction.atomic()
def recombine(genome, cassette, min_homology_arm_length,
              genome_name=None, cassette_name=None, notes=None):
    cassette = str(Seq(cassette))  # clean the sequence

    regions = find_swap_region(genome, cassette, min_homology_arm_length)
    if regions is None:
        return None

    new_genome = genome.update()
    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s recombined with %d bps integration cassette" % (genome.name,
                                                                          len(cassette))
    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    params = dict(cassette=cassette, homology_arm_length=min_homology_arm_length)
    op = Operation(genome=new_genome, type=Operation.RECOMBINATION[0], params=json.dumps(params))
    op.save()

    cassette_name = "Recombination cassette" if cassette_name is None else cassette_name

    # sort regions by starting positions, then by length. sorting by starting
    # position is important for shift_regions to work. sorting by length allows
    # shift_region to remove longer of the overlapping regions.
    regions = sorted(regions, key=lambda r: (r.start, len(r.sequence)))

    need_new_fragment = True
    while len(regions) > 0:
        region = regions[0]
        start, replaced, added, new_fragment_id =\
            recombine_region(new_genome, region, cassette,
                             min_homology_arm_length, cassette_name, op,
                             need_new_fragment)
        need_new_fragment = False
        regions = shift_regions(regions[1:],
                                region.fragment_id, start, replaced, added, new_fragment_id)

    return new_genome
