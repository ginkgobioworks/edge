import json
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
        pam_start = blast_res.subject_end + 1
        pam_end = pam_start + len(pam) - 1
        query = fragment.get_sequence(bp_lo=pam_start, bp_hi=pam_end)

    else:
        pam_end = blast_res.subject_end - 1
        pam_start = pam_end - len(pam) + 1
        query = fragment.get_sequence(bp_lo=pam_start, bp_hi=pam_end)
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


def crispr_dsb(genome, guide, pam, genome_name=None, notes=None):

    targets = find_crispr_target(genome, guide, pam)

    if len(targets) == 0:
        return None

    if genome_name is None or genome_name.strip() == "":
        genome_name = "%s after CRISPR-Cas9 wt (double stranded break) using guide %s"\
                      % (genome.name, guide)

    new_genome = genome.update()
    new_genome.name = genome_name
    new_genome.notes = notes
    new_genome.save()

    op = CrisprOp.get_operation(guide=guide, pam=pam)
    op.genome = new_genome
    op.save()

    for target in targets:
        if target.strand > 0:
            annotation_start = target.subject_start
            annotation_end = target.subject_end
        else:
            annotation_start = target.subject_end
            annotation_end = target.subject_start

        new_fragment_id = None
        with new_genome.update_fragment_by_fragment_id(target.fragment_id) as f:
            new_fragment_id = f.id

        with new_genome.annotate_fragment_by_fragment_id(new_fragment_id) as f:
            feature = 'CRISPR-Cas9 (pam %s) target' % pam
            f.annotate(annotation_start, annotation_end, feature,
                       'event', target.strand, operation=op)

    return new_genome


class CrisprOp(object):

    @staticmethod
    def check(genome, guide, pam, genome_name=None, notes=None):
        return find_crispr_target(genome, guide, pam)

    @staticmethod
    def get_operation(guide, pam, genome_name=None, notes=None):
        params = dict(guide=guide, pam=pam)
        op = Operation(type=Operation.CRISPR_DSB[0], params=json.dumps(params))
        return op

    @staticmethod
    def perform(genome, guide, pam, genome_name, notes):
        return crispr_dsb(genome, guide, pam, genome_name=genome_name, notes=notes)
