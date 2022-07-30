from edge.ssr import Reaction, Integration, Excision, Inversion, RMCE


class Sites(object):
    Frt  = _c("GAAGTTCCTATTC tctagaaa GTATAGGAACTTC")
    Frt1 = _c("GAAGTTCCTATTC tctagata GTATAGGAACTTC")
    Frt2 = _c("GAAGTTCCTATTC tctactta GTATAGGAACTTC")
    Frt3 = _c("GAAGTTCCTATTC ttcaaata GTATAGGAACTTC")
    Frt4 = _c("GAAGTTCCTATTC tctagaag GTATAGGAACTTC")
    Frt5 = _c("GAAGTTCCTATTC ttcaaaag GTATAGGAACTTC")


class FlpReaction(Reaction):

    @staticmethod
    def allowed():
        return [
        ]
