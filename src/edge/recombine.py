import os
import re
import json

from Bio.Seq import Seq

from edge.blast import blast_genome
from edge.models import Fragment, Operation
from edge.primer import design_primers_from_template
from edge.pcr import pcr_from_genome
from edge.orfs import detect_orfs


SINGLE_CROSSOVER_MAX_GAP = 2000
MAX_INTEGRATION_SIZE = 50000
CHECK_JUNCTION_PRIMER_WINDOW = 200
END_BPS_IGNORE = 10

CHECK_JUNCTION_LEFT_UP = 10
CHECK_JUNCTION_LEFT_DN = 190
CHECK_JUNCTION_RIGHT_UP = 190
CHECK_JUNCTION_RIGHT_DN = 10


def remove_overhangs(s):
    if s is None or len(s) == 0:
        return s
    s = re.sub(r"^\([\w\-\/]*\)", "", s)
    s = re.sub(r"\([\w\-\/]*\)$", "", s)
    s = re.sub(r"^\<[\w\-\/]*\>", "", s)
    s = re.sub(r"\<[\w\-\/]*\>$", "", s)
    return s


def get_cassette_new_annotations(
    cassette, cassette_reversed, specified_cassette, specified_annotations
):
    """
    user specified cassette sequence is `specified_cassette`, with
    `specified_annotations`. new sequence integrated into the genome is
    `cassette`. here we have to translate `specified_annotations` to
    coordinates relative to `cassette`.
    """

    cassette = cassette.lower()
    specified_cassette = specified_cassette.lower()

    annotations = []
    for orf in detect_orfs(cassette):
        annotations.append(
            dict(
                base_first=orf["start"],
                base_last=orf["end"],
                feature_name=orf["name"],
                feature_type="ORF",
                feature_strand=orf["strand"],
                feature_qualifiers=None,
            )
        )

    if specified_annotations:
        if cassette_reversed is True:
            cassette_len = len(cassette)
            reversed_cassette = str(Seq(cassette).reverse_complement())
            if reversed_cassette in specified_cassette:
                offset = specified_cassette.index(reversed_cassette)
                for annotation in specified_annotations:
                    if annotation["base_first"] - offset > 0:
                        annotations.append(
                            dict(
                                base_first=cassette_len
                                - (annotation["base_last"] - offset)
                                + 1,
                                base_last=cassette_len
                                - (annotation["base_first"] - offset)
                                + 1,
                                feature_name=annotation["name"],
                                feature_type=annotation["type"],
                                feature_strand=-annotation["strand"],
                                feature_qualifiers=annotation.get("qualifiers"),
                            )
                        )

        else:
            if cassette in specified_cassette:
                offset = specified_cassette.index(cassette)
                for annotation in specified_annotations:
                    if annotation["base_first"] - offset > 0 or True:
                        annotations.append(
                            dict(
                                base_first=annotation["base_first"] - offset,
                                base_last=annotation["base_last"] - offset,
                                feature_name=annotation["name"],
                                feature_type=annotation["type"],
                                feature_strand=annotation["strand"],
                                feature_qualifiers=annotation.get("qualifiers"),
                            )
                        )

    return annotations


class RecombinationRegion(object):
    def __init__(
        self,
        fragment_id,
        fragment_name,
        start,
        end,
        sequence,
        is_double_crossover,
        integrated_cassette,
        cassette_reversed,
        front_arm,
        back_arm,
    ):

        self.fragment_id = fragment_id
        self.fragment_name = fragment_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.is_double_crossover = is_double_crossover
        self.cassette_reversed = cassette_reversed
        if self.cassette_reversed:
            self.cassette = str(Seq(integrated_cassette).reverse_complement())
        else:
            self.cassette = integrated_cassette
        self.front_arm = front_arm
        self.back_arm = back_arm
        self.verification_front = []
        self.verification_back = []
        self.verification_cassette = []

    def to_dict(self):
        return self.__dict__


def compute_swap_region_from_results(
    front_arm_sequence,
    front_arm_blastres,
    back_arm_sequence,
    back_arm_blastres,
    cassette,
):
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
    if (
        front_arm_blastres.identity_ratio() < MIN_IDENTITIES
        or back_arm_blastres.identity_ratio() < MIN_IDENTITIES
        or front_arm_blastres.alignment_length()
        < len(front_arm_sequence) - MAX_MISSING_BP
        or back_arm_blastres.alignment_length()
        < len(back_arm_sequence) - MAX_MISSING_BP
    ):
        return None

    fragment = front_arm_blastres.fragment.indexed_fragment()

    single_crossover = False

    front_0 = front_arm_blastres.subject_start
    front_1 = front_arm_blastres.subject_end
    back_0 = back_arm_blastres.subject_start
    back_1 = back_arm_blastres.subject_end

    if front_arm_blastres.query_start > 1:
        cassette = cassette[front_arm_blastres.query_start - 1 :]
    if back_arm_blastres.query_end < len(back_arm_sequence):
        cassette = cassette[: -(len(back_arm_sequence) - back_arm_blastres.query_end)]

    # must align in the right orientation
    if front_arm_blastres.strand() > 0:  # on +1 strand
        is_reversed = False

        if (
            fragment.circular is False
            and back_1 < front_0
            and front_0 - back_1 <= SINGLE_CROSSOVER_MAX_GAP
        ) or (
            fragment.circular is True
            and (
                (back_1 < front_0 and front_0 - back_1 <= SINGLE_CROSSOVER_MAX_GAP)
                or (
                    back_1 < front_0 + fragment.length
                    and front_0 + fragment.length - back_1 <= SINGLE_CROSSOVER_MAX_GAP
                )
            )
        ):
            single_crossover = True
            region_start = back_0
            region_end = front_1
            single_crossover_front_dup_1 = front_0 - 1
            single_crossover_back_dup_0 = back_1 + 1

        elif (
            fragment.circular is False
            and front_1 < back_0
            and back_0 - front_1 <= MAX_INTEGRATION_SIZE
        ) or (
            fragment.circular is True
            and (
                (front_1 < back_0 and back_0 - front_1 <= MAX_INTEGRATION_SIZE)
                or (
                    front_1 < back_0 + fragment.length
                    and back_0 + fragment.length - front_1 <= MAX_INTEGRATION_SIZE
                )
            )
        ):
            single_crossover = False
            region_start = front_0
            region_end = back_1

        else:
            return None

    else:  # on -1 strand
        is_reversed = True

        if (
            fragment.circular is False
            and back_1 > front_0
            and back_1 - front_0 <= SINGLE_CROSSOVER_MAX_GAP
        ) or (
            fragment.circular is True
            and (
                (back_1 > front_0 and back_1 - front_0 <= SINGLE_CROSSOVER_MAX_GAP)
                or (
                    back_1 + fragment.length > front_0
                    and back_1 + fragment.length - front_0 <= SINGLE_CROSSOVER_MAX_GAP
                )
            )
        ):
            single_crossover = True
            region_start = front_1
            region_end = back_0
            single_crossover_front_dup_1 = back_1 - 1
            single_crossover_back_dup_0 = front_0 + 1

        elif (
            fragment.circular is False
            and front_1 > back_0
            and front_1 - back_0 <= MAX_INTEGRATION_SIZE
        ) or (
            fragment.circular is True
            and (
                (front_1 > back_0 and front_1 - back_0 <= MAX_INTEGRATION_SIZE)
                or (
                    front_1 + fragment.length > back_0
                    and front_1 + fragment.length - back_0 <= MAX_INTEGRATION_SIZE
                )
            )
        ):
            single_crossover = False
            region_start = back_1
            region_end = front_0

        else:
            return None

    if single_crossover is True:
        region_start = fragment.circ_bp(region_start)
        region_end = fragment.circ_bp(region_end)
        region = fragment.get_sequence(bp_lo=region_start, bp_hi=region_end)
        single_crossover_front_dup_1 = fragment.circ_bp(single_crossover_front_dup_1)
        single_crossover_back_dup_0 = fragment.circ_bp(single_crossover_back_dup_0)
        front_dup = fragment.get_sequence(
            bp_lo=region_start, bp_hi=single_crossover_front_dup_1
        )
        back_dup = fragment.get_sequence(
            bp_lo=single_crossover_back_dup_0, bp_hi=region_end
        )

        if is_reversed:
            cassette = (
                str(Seq(back_dup).reverse_complement())
                + cassette
                + str(Seq(front_dup).reverse_complement())
            )

        else:
            cassette = front_dup + cassette + back_dup

    else:
        # get sequence between arms, including arm
        region_start = fragment.circ_bp(region_start)
        region_end = fragment.circ_bp(region_end)
        region = fragment.get_sequence(bp_lo=region_start, bp_hi=region_end)

    is_double_crossover = not single_crossover
    return RecombinationRegion(
        fragment.id,
        fragment.name,
        region_start,
        region_end,
        region,
        is_double_crossover,
        cassette,
        is_reversed,
        front_arm_sequence,
        back_arm_sequence,
    )


def find_possible_swap_regions(genome, front_arm_sequence, back_arm_sequence, cassette):
    """
    Computes all possible recombined region from arm sequences.
    """

    front_arm_results = blast_genome(genome, "blastn", front_arm_sequence)
    back_arm_results = blast_genome(genome, "blastn", back_arm_sequence)

    regions = []
    for a_res in front_arm_results:
        for b_res in back_arm_results:
            region = compute_swap_region_from_results(
                front_arm_sequence, a_res, back_arm_sequence, b_res, cassette
            )
            if region is not None:
                regions.append(region)

    return regions


def remove_working_primers(genome, primers):
    failed_primers = []
    for primer in primers:
        p1 = primer["PRIMER_LEFT_SEQUENCE"]
        p2 = primer["PRIMER_RIGHT_SEQUENCE"]
        p = pcr_from_genome(genome, p1, p2)
        if p[0] is None:
            failed_primers.append(primer)
    return failed_primers


def get_verification_primers(genome, region, primer3_opts):
    """
    Design primers to verify replacing region with cassette.
    """

    fragment = Fragment.objects.get(pk=region.fragment_id).indexed_fragment()
    cassette = region.cassette

    check_junction_lu = CHECK_JUNCTION_LEFT_UP
    check_junction_ld = CHECK_JUNCTION_LEFT_DN
    check_junction_ru = CHECK_JUNCTION_RIGHT_UP
    check_junction_rd = CHECK_JUNCTION_RIGHT_DN

    if primer3_opts:
        if "EDGE_CHECK_JUNCTION_LEFT_HA" in primer3_opts:
            check_junction_ld = int(primer3_opts["EDGE_CHECK_JUNCTION_LEFT_HA"])
        if "EDGE_CHECK_JUNCTION_RIGHT_HA" in primer3_opts:
            check_junction_ru = int(primer3_opts["EDGE_CHECK_JUNCTION_RIGHT_HA"])

    #
    # front junction
    #

    upstream_window = CHECK_JUNCTION_PRIMER_WINDOW + check_junction_lu
    downstream_window = min(
        len(cassette), CHECK_JUNCTION_PRIMER_WINDOW + check_junction_ld
    )
    back = cassette[0:downstream_window]

    # don't bother looking for primers if front part of cassette is same as
    # region to be replaced
    if (
        region.start > 1
        and back.lower() != region.sequence[0:downstream_window].lower()
    ):
        front = fragment.get_sequence(region.start - upstream_window, region.start - 1)
        template = front + back
        junction = [len(front)]
        roi_start = (
            len(front) - check_junction_lu if len(front) > check_junction_lu else 0
        )
        roi_len = min(
            check_junction_lu + check_junction_ld,
            min(len(front), check_junction_ld) + len(cassette),
        )
        # remove primers that works on un-modified genome for creating a PCR product
        region.verification_front = remove_working_primers(
            genome,
            design_primers_from_template(
                template, roi_start, roi_len, junction, primer3_opts
            ),
        )

    #
    # back junction
    #

    upstream_window = min(
        len(cassette), CHECK_JUNCTION_PRIMER_WINDOW + check_junction_ru
    )
    downstream_window = CHECK_JUNCTION_PRIMER_WINDOW + check_junction_rd
    front = cassette[-upstream_window:]

    # don't bother looking for primers if front part of cassette is same as
    # region to be replaced
    if (
        region.start + len(region.sequence) < fragment.length
        and front.lower() != region.sequence[-upstream_window:].lower()
    ):
        back = fragment.get_sequence(
            region.start + len(region.sequence),
            region.start + len(region.sequence) + downstream_window - 1,
        )
        template = front + back
        junction = [len(front)]
        roi_start = (
            len(front) - check_junction_ru if len(front) > check_junction_ru else 0
        )
        roi_len = min(
            check_junction_ru + check_junction_rd,
            min(len(back), check_junction_rd) + len(cassette),
        )
        # remove primers that works on un-modified genome for creating a PCR product
        region.verification_back = remove_working_primers(
            genome,
            design_primers_from_template(
                template, roi_start, roi_len, junction, primer3_opts
            ),
        )

    #
    # cassette
    #

    if region.start > 1 and region.start + len(region.sequence) < fragment.length:
        front = fragment.get_sequence(
            region.start - CHECK_JUNCTION_PRIMER_WINDOW, region.start - 1
        )
        back = fragment.get_sequence(
            region.start + len(region.sequence),
            region.start + len(region.sequence) + CHECK_JUNCTION_PRIMER_WINDOW - 1,
        )
        template = front + cassette + back
        junctions = [len(front), len(front + cassette) - 1]
        roi_start = template.index(cassette)
        roi_len = len(cassette)
        region.verification_cassette = design_primers_from_template(
            template, roi_start, roi_len, junctions, primer3_opts
        )


def find_swap_region(
    genome,
    cassette,
    min_homology_arm_length,
    design_primers=False,
    primer3_opts=None,
    try_smaller_sequence=True,
):
    """
    Find a region on genome that can be recombined out using the cassette.
    Returns homology arms and possible regions along with cassette used to find matches.
    """

    cassette = remove_overhangs(cassette)
    if len(cassette) < 2 * min_homology_arm_length:
        return None

    front_arm = cassette[0:min_homology_arm_length]
    back_arm = cassette[-min_homology_arm_length:]

    regions = find_possible_swap_regions(genome, front_arm, back_arm, cassette)

    if len(regions) == 0 and try_smaller_sequence:
        cassette = cassette[END_BPS_IGNORE:-END_BPS_IGNORE]
        if len(cassette) >= 2 * min_homology_arm_length:
            return find_swap_region(
                genome,
                cassette,
                min_homology_arm_length,
                design_primers=design_primers,
                primer3_opts=primer3_opts,
                try_smaller_sequence=False,
            )
        return []

    if design_primers is True:
        for region in regions:
            get_verification_primers(genome, region, primer3_opts)

    return regions


def recombine_region(genome, region, min_homology_arm_length, op, new_fragment_dict):
    """
    Recombines on a given region. Returns recombination cassette location, how
    many base pairs removed, how many base pairs added, new fragment id.
    new_fragment_dict is a dictionary that indicates if the region.fragment_id is already
    a newly created fragment.
    """

    region_start = region.start
    region_end = region.end
    cassette = region.cassette

    new_fragment_id = None
    is_new_fragment = region.fragment_id not in new_fragment_dict

    with genome.update_fragment_by_fragment_id(
        region.fragment_id, new_fragment=is_new_fragment
    ) as f:
        new_fragment_dict[f.id] = 1

        replaced = 0

        # added "+1 here": region_start and region_end are interpreted as
        # replacing the bp starting at region_start and ending at region_end.
        # If region_end is region_start-1, that means it's just an insertion
        # and we should allow that.

        if region_start <= region_end + 1:
            number_of_bps_to_remove = region_end - region_start + 1
            f.replace_bases(region_start, number_of_bps_to_remove, cassette)
            replaced = region_end - region_start + 1
        else:
            assert f.circular is True
            f.replace_bases(region_start, f.length - region_start + 1, cassette)
            replaced = f.length - region_start + 1
            f.remove_bases(1, region_end)
            replaced += region_end
            # adjust region_start after removing sequence at the start
            region_start -= region_end

        new_fragment_id = f.id

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

    for r in regions:
        if r.fragment_id == fragment_id:
            assert start <= r.start

    def updated_region(region):
        if region.start <= start + replaced - 1:  # overlapping
            return None
        region.fragment_id = new_fragment_id
        if region.start < region.end:
            region.start = region.start - replaced + added
            region.end = region.end - replaced + added
        else:
            # last replacement cannot be across circular boundary, since we are
            # replacing in order of starting position
            region.start = region.start - replaced + added
        return region

    regions = [
        region if region.fragment_id != fragment_id else updated_region(region)
        for region in regions
    ]
    return [r for r in regions if r is not None]


def find_common_start(a, b):
    return len(os.path.commonprefix([a, b]))


def trim_cassette_and_region(genome, region):
    """
    Trim an integration region by removing sequences at both ends that also
    appear on the specified integration donor cassette (i.e. removing homology
    arms).
    """

    if region.is_double_crossover is False:
        return region

    fragment = [f for f in genome.fragments.all() if f.id == region.fragment_id][0]
    fragment_length = fragment.indexed_fragment().length

    integrated_sequence = region.cassette
    integrated_sequence = integrated_sequence.lower()
    replaced_sequence = region.sequence.lower()
    if integrated_sequence == replaced_sequence:
        return None

    # print("region %s-%s, fragment length %s" % (region.start, region.end, fragment_length))
    # print("integrated: %s" % integrated_sequence)
    # print("  original: %s" % replaced_sequence)

    trim_front = find_common_start(integrated_sequence, replaced_sequence)
    if trim_front:
        region.cassette = region.cassette[trim_front:]
        # print("trim front", trim_front)
        region.start += trim_front
        region.start = ((region.start - 1) % fragment_length) + 1

        integrated_sequence = integrated_sequence[trim_front:]
        replaced_sequence = replaced_sequence[trim_front:]

    trim_back = find_common_start(integrated_sequence[::-1], replaced_sequence[::-1])
    if trim_back:
        region.cassette = region.cassette[:-trim_back]
        # print("trim back", trim_back)
        region.end -= trim_back
        region.end = ((region.end - 1) % fragment_length) + 1

    # print("new coordinates %s-%s" % (region.start, region.end))
    # print("new cassette %s" % region.cassette)
    return region


def recombine_sequence(
    genome,
    cassette,
    homology_arm_length,
    genome_name=None,
    cassette_name=None,
    notes=None,
    annotations=None
):
    cassette = remove_overhangs(cassette)
    cassette = str(Seq(cassette))  # clean the sequence

    regions = find_swap_region(genome, cassette, homology_arm_length)
    if regions is None or len(regions) == 0:
        return None

    regions = [trim_cassette_and_region(genome, r) for r in regions]
    regions = [r for r in regions if r is not None]
    if len(regions) == 0:
        return None

    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s recombined with %d bps integration cassette" % (
            genome.name,
            len(cassette),
        )
    new_genome = genome.update()
    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    op = RecombineOp.get_operation(
        cassette=cassette, homology_arm_length=homology_arm_length
    )
    op.genome = new_genome
    op.save()

    cassette_name = "Integrated sequence" if cassette_name is None else cassette_name

    # sort regions by starting positions, then by length. sorting by starting
    # position is important for shift_regions to work. sorting by length allows
    # shift_region to remove longer of the overlapping regions.
    regions = sorted(regions, key=lambda r: (r.start, len(r.sequence)))

    regions_before = [
        dict(
            fragment_id=region.fragment_id,
            start=region.start,
            end=region.end,
            cassette_reversed=region.cassette_reversed,
            cassette=region.cassette
        )
        for region in regions
    ]
    regions_after = []

    new_fragment_dict = {}

    while len(regions) > 0:
        region = regions[0]

        # passing new_fragment_dict to recombine_region, so we can remember which
        # fragment is already newly created
        start, replaced, added, new_fragment_id = recombine_region(
            new_genome, region, homology_arm_length, op, new_fragment_dict
        )

        regions_after.append(dict(fragment_id=region.fragment_id, start=start))

        # after first integration, we need to adjust computed coordinates for
        # other loci by effects of the first integration
        regions = shift_regions(
            regions[1:], region.fragment_id, start, replaced, added, new_fragment_id
        )
        for r in regions_after:
            if r["fragment_id"] == region.fragment_id:
                r["fragment_id"] = new_fragment_id

    return dict(
        new_genome=new_genome,
        regions=dict(before=regions_before, after=regions_after),
        cassette_name=cassette_name,
        operation=op,
    )


def annotate_integration(
    genome,
    new_genome,
    regions_before,
    regions_after,
    cassette_name,
    op,
    specified_cassette,
    specified_annotations,
):

    specified_cassette = specified_cassette.lower()
    reversed_specified_cassette = str(
        Seq(specified_cassette).reverse_complement()
    ).lower()

    before_and_after_with_annotations = []
    for before, after in zip(regions_before, regions_after):
        annotations_to_use = []
        if (
            before["cassette"].lower() not in specified_cassette
            and before["cassette"].lower() not in reversed_specified_cassette
        ):
            # opps, annotations no longer match, ignore them for now
            annotations_to_use = []
        else:
            annotations_to_use = specified_annotations

        annotations = get_cassette_new_annotations(
            before["cassette"],
            before["cassette_reversed"],
            specified_cassette,
            annotations_to_use,
        )
        before_and_after_with_annotations.append((before, after, annotations))

    for before, after, annotations in before_and_after_with_annotations:
        with new_genome.annotate_fragment_by_fragment_id(after["fragment_id"]) as f:
            # region_start is already adjusted for multiple integration in this
            # fragment

            # collect list of annotations to be added later

            # annotate cassette
            if before["cassette_reversed"]:
                strand = -1
            else:
                strand = 1
            if len(before["cassette"]) > 0:
                f.annotate(
                    after["start"],
                    after["start"] + len(before["cassette"]) - 1,
                    cassette_name,
                    "operation",
                    strand,
                    operation=op,
                )

            # annotated inherited and new annotations
            for annotation in annotations:
                f.annotate(
                    after["start"] + annotation["base_first"] - 1,
                    after["start"] + annotation["base_last"] - 1,
                    annotation["feature_name"],
                    annotation["feature_type"],
                    annotation["feature_strand"],
                    qualifiers=annotation.get("feature_qualifiers"),
                )


def recombine(
    genome,
    cassette,
    homology_arm_length,
    genome_name=None,
    cassette_name=None,
    notes=None,
    annotations=None,
):

    # lock root genome to prevent other genomes of touching same fragment or chunk
    genome.lock()

    x = recombine_sequence(
        genome,
        cassette,
        homology_arm_length,
        genome_name=genome_name,
        cassette_name=cassette_name,
        notes=notes,
        annotations=annotations
    )

    if x is None:
        return x
    new_genome = x["new_genome"]

    annotate_integration(
        genome,
        new_genome,
        x["regions"]["before"],
        x["regions"]["after"],
        x["cassette_name"],
        x["operation"],
        cassette,
        annotations
    )

    return new_genome


class RecombineOp(object):
    @staticmethod
    def check(
        genome,
        cassette,
        homology_arm_length,
        genome_name=None,
        cassette_name=None,
        notes=None,
        design_primers=False,
        primer3_opts=None,
        annotations=None,
    ):
        return find_swap_region(
            genome,
            cassette,
            homology_arm_length,
            design_primers=design_primers,
            primer3_opts=primer3_opts,
        )

    @staticmethod
    def get_operation(
        cassette,
        homology_arm_length,
        genome_name=None,
        cassette_name=None,
        notes=None,
        design_primers=False,
        primer3_opts=None,
        annotations=None,
    ):
        cassette = remove_overhangs(cassette)
        params = dict(cassette=cassette, homology_arm_length=homology_arm_length)
        op = Operation(type=Operation.RECOMBINATION[0], params=json.dumps(params))
        return op

    @staticmethod
    def perform(
        genome,
        cassette,
        homology_arm_length,
        genome_name,
        cassette_name,
        notes,
        design_primers=False,
        primer3_opts=None,
        annotations=None,
    ):
        new_genome = recombine(
            genome,
            cassette,
            homology_arm_length,
            genome_name=genome_name,
            cassette_name=cassette_name,
            notes=notes,
            annotations=annotations,
        )

        # some fragment indices may be destroyed because we had to split up
        # existing chunks; rebuild indices so that parent genome is immediately
        # reusable in an RO server
        genome.indexed_genome()

        return new_genome
