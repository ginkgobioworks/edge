from edge.ssr import rc, _c, Reaction, Integration, Excision, Inversion, RMCE


# TODO:
#
# See original SSR Edge dev plan
# https://docs.google.com/document/d/15iXSn6uWsn16v7nkkxzP9r2korEwr39G_sRlYKbX6IY/edit?pli=1
#
# 1. Add attL and attR sites
# 2. Add rules for integration and RMCE, allowing attB_rc(attP) as well as attB_attP
# 3. Add alternative/orthogonal attB/P sites
# 4. Add separate PhiC31 and Bxb1 reactions
# 5. Figure out how often they would be used in combo, and if we need combined reaction


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
