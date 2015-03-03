import json
from django.db import transaction
from edge.blast import blast_genome
from edge.models import Fragment, Operation
from edge.primer import design_primers_from_template
from edge.pcr import pcr_from_genome
from Bio.Seq import Seq


CHECK_JUNCTION_PRIMER_WINDOW = 200
CHECK_JUNCTION_SIZE = 200
END_BPS_IGNORE = 8


def remove_overhangs(s):
    if s is None or len(s) == 0:
        return s
    if s[0] == '(' and s.find(')') >= 0:
        s = s[s.find(')')+1:]
    if s[-1] == ')' and s.rfind('(') >= 0:
        s = s[:s.rfind('(')]
    return s


def blast_result_annotations(res):
    fragment = Fragment.objects.get(pk=res.fragment_id).indexed_fragment()

    # get annotations, but only those that correspond to a full feature
   
    annotations = [dict(base_first=a.base_first,
                        base_last=a.base_last,
                        feature_name=a.feature.name,
                        fragment_id=fragment.id)
                   for a in fragment.annotations(res.subject_start, res.subject_end)
                   if a.feature_base_first == 1 and a.feature_base_last == a.feature.length and
                      a.base_first >= res.subject_start and a.base_last <= res.subject_end]

    return annotations


def identify_changes(query, subject):
    # identify changes w/r/t subject, assuming subject is complete sequence of
    # some feature. query and subject are alignments from blast.

    diffs = []
    si = 0

    for i, s in enumerate(subject):
        s = s.upper()
        q = query[i].upper()

        if s == '-':
            # insertion in query
            diffs.append('+%s%s' % (si, q))

        elif q == '-':
            # deletion in query
            diffs.append('-%s%s' % (si, s))
            si += 1

        elif q != s:
            # mutation
            diffs.append('%s%s%s' % (s, si, q))
            si += 1

        else:
            # same
            si += 1

    return diffs


def get_cassette_annotation(annotation, blast_res):
    # for an annotation in blast result, returns fragment of query from blast
    # result and differences to subject

    q = []
    s = []

    subject_i = blast_res.subject_start-1
    query_i = blast_res.query_start-1
    query_start = None
    query_end = None

    for i, subject_base in enumerate(blast_res.alignment['subject']):
        if subject_base != '-':
            subject_i += 1
      
        query_base = blast_res.alignment['query'][i]
        if query_base != '-':
            query_i += 1

        if subject_i >= annotation['base_first'] and subject_i <= annotation['base_last']:
            if subject_i == annotation['base_first']:
                query_start = query_i
            if subject_i == annotation['base_last']:
                query_end = query_i
            q.append(blast_res.alignment['query'][i])
            s.append(blast_res.alignment['subject'][i])

    diffs = identify_changes(q, s)
    diffs = ' '.join(diffs)

    return dict(base_first=query_start, base_last=query_end,
                feature_name=annotation['feature_name']+' '+diffs)


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
        self.verification_front = []
        self.verification_back = []
        self.verification_cassette = []

    def to_dict(self):
        return self.__dict__

    def update_annotations(self, genome, cassette):
        self.cassette_annotations = self.get_cassette_annotations(genome, cassette)

    def get_cassette_annotations(self, genome, cassette):
        return self.get_cassette_inherited_annotations(genome, cassette) +\
               self.get_cassette_new_annotations(genome, cassette)

    def get_cassette_inherited_annotations(self, genome, cassette):
        # blast cassette against unmodified genome
        if self.cassette_reversed:
            cassette = str(Seq(cassette).reverse_complement())
        cassette_blast_res = blast_genome(genome, 'blastn', cassette)

        # find matches inside region to be replaced
        matches = [res for res in cassette_blast_res\
                   if res.fragment_id == self.fragment_id and
                      res.subject_start >= self.start and res.subject_end <= self.end]

        inherited_annotations = []

        for blast_res in matches:
            # get annotations for each match
            annotations = blast_result_annotations(blast_res)

            print blast_res.to_dict()
            print annotations

            # update annotations with blast result
            for annotation in annotations:
                inherited_annotations.append(get_cassette_annotation(annotation, blast_res))

        return inherited_annotations

    def get_cassette_new_annotations(self, genome, cassette):
        # XXX
        # detect ORF on cassette

        return []


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


def remove_working_primers(genome, primers):
    failed_primers = []
    for primer in primers:
        p1 = primer['PRIMER_LEFT_SEQUENCE']
        p2 = primer['PRIMER_RIGHT_SEQUENCE']
        p = pcr_from_genome(genome, p1, p2)
        if p[0] is None:
            failed_primers.append(primer)
    return failed_primers


def get_verification_primers(genome, region, cassette, primer3_opts):
    """
    Design primers to verify replacing region with cassette.
    """

    fragment = Fragment.objects.get(pk=region.fragment_id).indexed_fragment()
    if region.cassette_reversed:
        cassette = str(Seq(cassette).reverse_complement())

    #
    # front junction
    #

    upstream_window = CHECK_JUNCTION_PRIMER_WINDOW+CHECK_JUNCTION_SIZE/2
    downstream_window = min(len(cassette), CHECK_JUNCTION_PRIMER_WINDOW+CHECK_JUNCTION_SIZE/2)
    back = cassette[0:downstream_window]

    # don't bother looking for primers if front part of cassette is same as
    # region to be replaced
    if region.start > 1 and back.lower() != region.sequence[0:downstream_window].lower():
        front = fragment.get_sequence(region.start-upstream_window, region.start-1)
        template = front+back
        junction = [len(front)]
        roi_start = len(front)-CHECK_JUNCTION_SIZE/2 if len(front) > CHECK_JUNCTION_SIZE/2 else 0
        roi_len = min(CHECK_JUNCTION_SIZE, min(len(front), CHECK_JUNCTION_SIZE/2)+len(cassette))
        # remove primers that works on un-modified genome for creating a PCR product
        region.verification_front =\
            remove_working_primers(genome,
                                   design_primers_from_template(template, roi_start, roi_len,
                                                                junction, primer3_opts))

    #
    # back junction
    #

    upstream_window = min(len(cassette), CHECK_JUNCTION_PRIMER_WINDOW+CHECK_JUNCTION_SIZE/2)
    downstream_window = CHECK_JUNCTION_PRIMER_WINDOW+CHECK_JUNCTION_SIZE/2
    front = cassette[-upstream_window:]

    # don't bother looking for primers if front part of cassette is same as
    # region to be replaced
    if region.start+len(region.sequence) < fragment.length and\
       front.lower() != region.sequence[-upstream_window:].lower():
        back = fragment.get_sequence(region.start+len(region.sequence),
                                     region.start+len(region.sequence)+downstream_window-1)
        template = front+back
        junction = [len(front)]
        roi_start = len(front)-CHECK_JUNCTION_SIZE/2 if len(front) > CHECK_JUNCTION_SIZE/2 else 0
        roi_len = min(CHECK_JUNCTION_SIZE, min(len(back), CHECK_JUNCTION_SIZE/2)+len(cassette))
        # remove primers that works on un-modified genome for creating a PCR product
        region.verification_back =\
            remove_working_primers(genome,
                                   design_primers_from_template(template, roi_start, roi_len,
                                                                junction, primer3_opts))

    #
    # cassette
    #

    if region.start > 1 and region.start+len(region.sequence) < fragment.length:
        front = fragment.get_sequence(region.start-CHECK_JUNCTION_PRIMER_WINDOW, region.start-1)
        back = fragment.get_sequence(region.start+len(region.sequence),
                                     region.start+len(region.sequence) +
                                     CHECK_JUNCTION_PRIMER_WINDOW-1)
        template = front+cassette+back
        junctions = [len(front), len(front+cassette)-1]
        roi_start = template.index(cassette)
        roi_len = len(cassette)
        region.verification_cassette =\
            design_primers_from_template(template, roi_start, roi_len,
                                         junctions, primer3_opts)


def _find_swap_region(genome, cassette, min_homology_arm_length,
                      design_primers=False, primer3_opts=None, try_smaller_sequence=True):
    """
    Find a region on genome that can be recombined out using the cassette.
    Returns homology arms and possible regions along with cassette used to find matches.
    """

    cassette = remove_overhangs(cassette)
    if len(cassette) < 2*min_homology_arm_length:
        return (None, cassette)

    front_arm = cassette[0:min_homology_arm_length]
    back_arm = cassette[-min_homology_arm_length:]

    regions = find_possible_swap_regions(genome, front_arm, back_arm)
    if len(regions) == 0 and try_smaller_sequence:
        cassette = cassette[END_BPS_IGNORE:-END_BPS_IGNORE]
        if len(cassette) >= 2*min_homology_arm_length:
            return _find_swap_region(genome, cassette, min_homology_arm_length,
                                     design_primers, False)
        return ([], cassette)

    if design_primers is True:
        for region in regions:
            get_verification_primers(genome, region, cassette, primer3_opts)

    return (regions, cassette)


def find_swap_region(genome, cassette, min_homology_arm_length,
                     design_primers=False, primer3_opts=None):
    """
    Same as _find_swap_region, but does not return cassette used to find matches.
    """

    x = _find_swap_region(genome, cassette, min_homology_arm_length,
                          design_primers, primer3_opts, True)

    regions = []
    for region in x[0]:
        region.update_annotations(genome, cassette)
        regions.append(region)

    return regions


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
def recombine(genome, cassette, homology_arm_length,
              genome_name=None, cassette_name=None, notes=None):
    cassette = remove_overhangs(cassette)
    cassette = str(Seq(cassette))  # clean the sequence

    regions, cassette = _find_swap_region(genome, cassette, homology_arm_length)
    if regions is None:
        return None

    new_genome = genome.update()
    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s recombined with %d bps integration cassette" % (genome.name,
                                                                          len(cassette))
    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    op = RecombineOp.get_operation(cassette=cassette, homology_arm_length=homology_arm_length)
    op.genome = new_genome
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
                             homology_arm_length, cassette_name, op,
                             need_new_fragment)
        need_new_fragment = False
        regions = shift_regions(regions[1:],
                                region.fragment_id, start, replaced, added, new_fragment_id)

    return new_genome


class RecombineOp(object):

    @staticmethod
    def check(genome, cassette, homology_arm_length,
              genome_name=None, cassette_name=None, notes=None,
              design_primers=False, primer3_opts=None):
        return find_swap_region(genome, cassette, homology_arm_length,
                                design_primers=design_primers,
                                primer3_opts=primer3_opts)

    @staticmethod
    def get_operation(cassette, homology_arm_length,
                      genome_name=None, cassette_name=None, notes=None,
                      design_primers=False, primer3_opts=None):
        cassette = remove_overhangs(cassette)
        params = dict(cassette=cassette, homology_arm_length=homology_arm_length)
        op = Operation(type=Operation.RECOMBINATION[0], params=json.dumps(params))
        return op

    @staticmethod
    def perform(genome, cassette, homology_arm_length, genome_name, cassette_name, notes,
                design_primers=False, primer3_opts=None):
        return recombine(genome, cassette, homology_arm_length,
                         genome_name=genome_name, cassette_name=cassette_name, notes=notes)
