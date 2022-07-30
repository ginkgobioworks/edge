from edge.ssr import Reaction, Integration, Excision, Inversion, RMCE


class Sites(object):
    loxP     = _c("ATAACTTCGTATA GCATACAT TATACGAAGTTAT")
    lox66    = _c("ATAACTTCGTATA GCATACAT TATACGAACGGTA")
    lox71    = _c("TACCGTTCGTATA GCATACAT TATACGAAGTTAT")
    lox72    = _c("TACCGTTCGTATA GCATACAT TATACGAACGGTA")

    lox5171  = _c("ATAACTTCGTATA GtAcACAT TATACGAAGTTAT")
    lox2272  = _c("ATAACTTCGTATA GgATACtT TATACGAAGTTAT")

    loxm2    = _c("ATAACTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_66 = _c("ATAACTTCGTATA TGGTTTCT TATACGAACGGTA")
    loxm2_71 = _c("TACCGTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_72 = _c("TACCGTTCGTATA TGGTTTCT TATACGAACGGTA")


class CreLoxReaction(Reaction):

    @staticmethod
    def allowed():
        return [
            Integration(Sites.lox66, Sites.lox71, Sites.lox72, Sites.loxP),
            Integration(Sites.lox71, Sites.lox66, Sites.loxP, Sites.lox72),

            Excision(Sites.loxP, Sites.loxP, Sites.loxP),
            Excision(Sites.lox66, Sites.loxP, Sites.loxP),
            Excision(Sites.loxP, Sites.lox66, Sites.lox66),
            Excision(Sites.lox71, Sites.loxP, Sites.lox71),
            Excision(Sites.loxP, Sites.lox71, Sites.loxP),

            Inversion(Sites.loxP, rc(Sites.loxP), Sites.loxP, rc(Sites.loxP)),
            Inversion(Sites.lox66, rc(Sites.loxP), Sites.loxP, rc(Sites.lox66)),
            Inversion(Sites.loxP, rc(Sites.lox66), Sites.lox66, rc(Sites.loxP)),
            Inversion(Sites.lox71, rc(Sites.loxP), Sites.lox71, rc(Sites.loxP)),
            Inversion(Sites.loxP, rc(Sites.lox71), Sites.loxP, rc(Sites.lox71)),

            RMCE(Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272),
            RMCE(Sites.lox66, Sites.lox2272, Sites.lox71, Sites.lox2272, Sites.lox72, Sites.lox2272),
        ]
