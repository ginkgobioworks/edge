import json
from django.db import transaction
from edge.blast import blast_genome
from edge.models import Operation
from Bio.Seq import Seq


class CrisprTarget(object):

    def __init__(self, fragment_id, fragment_name, strand, subject_start, subject_end, pam):
        self.fragment_id = fragment_id
        self.fragment_name = fragment_name
        self.strand = strand
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.pam = pam

    def to_dict(self):
        return self.__dict__


def match_pam(pam, query):
    if len(query) != len(pam):
        return False
    for p, q in zip(pam.lower(), query.lower()):
        if p != 'n' and p != q:
            return False
    return True


def target_followed_by_pam(blast_res, pam):
    fragment = blast_res.fragment.indexed_fragment()

    if blast_res.strand() > 0:
        pam_start = blast_res.subject_end+1
        pam_end = pam_start+len(pam)-1

        pam_start = fragment.circ_bp(pam_start)
        pam_end = fragment.circ_bp(pam_end)

        if pam_start < pam_end:
            query = fragment.get_sequence(bp_lo=pam_start, bp_hi=pam_end)
        else:
            assert fragment.circular is True
            query_p1 = fragment.get_sequence(bp_lo=pam_start)
            query_p2 = fragment.get_sequence(bp_hi=pam_end)
            query = query_p1+query_p2

    else:
        pam_end = blast_res.subject_end-1
        pam_start = pam_end-len(pam)+1

        pam_end = fragment.circ_bp(pam_end)
        pam_start = fragment.circ_bp(pam_start)

        if pam_start < pam_end:
            query = fragment.get_sequence(bp_lo=pam_start, bp_hi=pam_end)
        else:
            assert fragment.circular is True
            query_p1 = fragment.get_sequence(bp_lo=pam_start)
            query_p2 = fragment.get_sequence(bp_hi=pam_end)
            query = query_p1+query_p2
        query = str(Seq(query).reverse_complement())

    if match_pam(pam, query) is True:
        subject_start = blast_res.subject_start
        subject_end = blast_res.subject_end
        subject_start = fragment.circ_bp(subject_start)
        subject_end = fragment.circ_bp(subject_end)

        return CrisprTarget(blast_res.fragment_id, blast_res.fragment.name,
                            blast_res.strand(), subject_start, subject_end, pam)
    return None


def find_crispr_target(genome, guide, pam):
    """
    Find sequences on genome that have exact match to guide, followed by pam
    sequence.
    """

    guide_matches = blast_genome(genome, 'blastn', guide)
    targets = []

    for res in guide_matches:
        if res.query_start == 1 and res.query_end == len(guide):
            target = target_followed_by_pam(res, pam)
            if target is not None:
                targets.append(target)

    return targets


@transaction.atomic()
def crispr_dsb(genome, guide, pam, cut_bps_before_pam, genome_name=None, notes=None):

    targets = find_crispr_target(genome, guide, pam)

    if len(targets) == 0:
        return None

    # we update 2 bases around the CRISPR cut site to Ns, to denote the
    # mutagenic nature of the repair
    bases_to_update = 2

    new_genome = genome.update()
    for target in targets:
        if target.strand > 0:
            repair_bp1 = target.subject_end-cut_bps_before_pam
            repair_bp2 = repair_bp1+bases_to_update-1
        else:
            repair_bp1 = target.subject_end+cut_bps_before_pam-1
            repair_bp2 = repair_bp1+bases_to_update-1

        fragment = genome.fragments.filter(id=target.fragment_id)[0].indexed_fragment()
        repair_bp1 = fragment.circ_bp(repair_bp1)
        repair_bp2 = fragment.circ_bp(repair_bp2)
        new_fragment_id = None

        with new_genome.update_fragment_by_fragment_id(target.fragment_id) as f:
            if repair_bp1 < repair_bp2:
                f.replace_bases(repair_bp1, bases_to_update, 'N'*bases_to_update)
            else:
                assert f.circular is True
                n = f.length-repair_bp1+1
                f.replace_bases(repair_bp1, n, 'N'*n)
                f.replace_bases(1, repair_bp2, 'N'*repair_bp2)

            new_fragment_id = f.id

        with new_genome.annotate_fragment_by_fragment_id(new_fragment_id) as f:
            if repair_bp1 < repair_bp2:
                f.annotate(repair_bp1, repair_bp2, 'CRISPR DSB repair', 'feature', target.strand)
            else:
                n = f.length-repair_bp1+1
                f.annotate(repair_bp1, n, 'CRISPR DSB repair', 'feature', target.strand)
                f.annotate(1, repair_bp2, 'CRISPR DSB repair', 'feature', target.strand)

    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s with CRISPR DSB repair, using guide %s" % (genome.name, guide)

    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    params = dict(guide=guide, pam=pam)
    op = Operation(type=Operation.CRISPR_DSB[0], params=json.dumps(params))
    op.save()
    new_genome.operations.add(op)

    return new_genome
