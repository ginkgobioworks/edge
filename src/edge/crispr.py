import json
from django.db import transaction
from edge.blast import blast_genome
from edge.models import Fragment, Operation
from Bio.Seq import Seq


class CrisprTarget(object):

    def __init__(self, fragment_id, fragment_name, target_start, target_end, pam):
        self.fragment_id = fragment_id
        self.fragment_name = fragment_name
        self.target_start = target_start
        self.target_end = target_end
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

        pam_start = ((pam_start-1)%fragment.length)+1
        pam_end = ((pam_end-1)%fragment.length)+1

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

        pam_end = ((pam_end-1)%fragment.length)+1
        pam_start = ((pam_start-1)%fragment.length)+1

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
        subject_start = ((subject_start-1)%fragment.length)+1
        subject_end = ((subject_end-1)%fragment.length)+1

        return CrisprTarget(blast_res.fragment_id, blast_res.fragment.name,
                            subject_start, subject_end, pam)
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
