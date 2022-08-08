import itertools
from edge.ssr import rc, lower_no_whitespace, Reaction, Excision, Inversion, RMCE


class Sites(object):
    Frt = lower_no_whitespace("GAAGTTCCTATTC tctagaaa GTATAGGAACTTC")
    Frt1 = lower_no_whitespace("GAAGTTCCTATTC tctagata GTATAGGAACTTC")
    Frt2 = lower_no_whitespace("GAAGTTCCTATTC tctactta GTATAGGAACTTC")
    Frt3 = lower_no_whitespace("GAAGTTCCTATTC ttcaaata GTATAGGAACTTC")
    Frt4 = lower_no_whitespace("GAAGTTCCTATTC tctagaag GTATAGGAACTTC")
    Frt5 = lower_no_whitespace("GAAGTTCCTATTC ttcaaaag GTATAGGAACTTC")


class FlpReaction(Reaction):

    @staticmethod
    def allowed():
        frt_sites = [Sites.Frt, Sites.Frt1, Sites.Frt2, Sites.Frt3, Sites.Frt4, Sites.Frt5]

        allowed = []
        for frt in frt_sites:
            allowed.append(Excision(frt, frt, frt))
            allowed.append(Inversion(frt, rc(frt), frt, rc(frt)))

        for combo in list(itertools.combinations(frt_sites, 2)):
            allowed.append(RMCE(combo[0], combo[1], combo[0], combo[1], combo[0], combo[1]))
            allowed.append(RMCE(combo[1], combo[0], combo[1], combo[0], combo[1], combo[0]))

        return allowed
