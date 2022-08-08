import itertools
from edge.ssr import rc, _c, Reaction, Excision, Inversion, RMCE


class Sites(object):
    Frt = _c("GAAGTTCCTATTC tctagaaa GTATAGGAACTTC")
    Frt1 = _c("GAAGTTCCTATTC tctagata GTATAGGAACTTC")
    Frt2 = _c("GAAGTTCCTATTC tctactta GTATAGGAACTTC")
    Frt3 = _c("GAAGTTCCTATTC ttcaaata GTATAGGAACTTC")
    Frt4 = _c("GAAGTTCCTATTC tctagaag GTATAGGAACTTC")
    Frt5 = _c("GAAGTTCCTATTC ttcaaaag GTATAGGAACTTC")


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
