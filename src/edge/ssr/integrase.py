from edge.ssr import rc, _c, Reaction, Integration, Excision, Inversion, RMCE


class Sites(object):
    """
    PhiC31 sites
    """

    PhiC31_attB_tt = _c("TGCGGGTGCCAGGGCGTGCCC tt GGGCTCCCCGGGCGCGTACTCC")
    PhiC31_attP_tt = _c("GTGCCCCAACTGGGGTAACCT tt GAGTTCTCTCAGTTGGGGG")

    """
    Bxb1 sites
    """

    Bxb1_attB =     _c("TCGGCCGGCTTGTCGACGACG gcggtctc CGTCGTCAGGATCATCCGGGC")
    Bxb1_attP = _c("GTCGTGGTTTGTCTGGTCAACCACC gcggtctc AGTGGTGTACGGTACAAACCCCGAC")


class Bxb1Reaction(Reaction):

    @staticmethod
    def allowed():
        return [
        ]


class PhiC31Reaction(Reaction):

    @staticmethod
    def allowed():
        return [
        ]
