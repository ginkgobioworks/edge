from edge.ssr import rc, lower_no_whitespace, Reaction, Integration, Excision, Inversion, RMCE


class Sites(object):
    loxP = lower_no_whitespace("ATAACTTCGTATA GCATACAT TATACGAAGTTAT")
    lox66 = lower_no_whitespace("ATAACTTCGTATA GCATACAT TATACGAACGGTA")
    lox71 = lower_no_whitespace("TACCGTTCGTATA GCATACAT TATACGAAGTTAT")
    lox72 = lower_no_whitespace("TACCGTTCGTATA GCATACAT TATACGAACGGTA")

    lox5171 = lower_no_whitespace("ATAACTTCGTATA GtAcACAT TATACGAAGTTAT")
    lox2272 = lower_no_whitespace("ATAACTTCGTATA GgATACtT TATACGAAGTTAT")

    loxm2 = lower_no_whitespace("ATAACTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_66 = lower_no_whitespace("ATAACTTCGTATA TGGTTTCT TATACGAACGGTA")
    loxm2_71 = lower_no_whitespace("TACCGTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_72 = lower_no_whitespace("TACCGTTCGTATA TGGTTTCT TATACGAACGGTA")


class CreLoxReaction(Reaction):

    @staticmethod
    def allowed():
        return [
            RMCE(Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272),
            RMCE(Sites.lox2272, Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272, Sites.loxP),
            RMCE(Sites.lox66, Sites.lox2272,
                 Sites.lox71, Sites.lox2272,
                 Sites.lox72, Sites.lox2272),
            RMCE(Sites.lox2272, Sites.lox66,
                 Sites.lox2272, Sites.lox71,
                 Sites.lox2272, Sites.lox72),
            RMCE(Sites.loxm2_71, Sites.lox66,
                 Sites.loxm2_66, Sites.lox71,
                 Sites.loxP, Sites.lox72),
            RMCE(Sites.loxm2_66, Sites.lox66,
                 Sites.loxm2_71, Sites.lox71,
                 Sites.loxm2_72, Sites.lox72),

            Integration(Sites.lox66, Sites.lox71, Sites.lox72, Sites.loxP),
            Integration(Sites.lox71, Sites.lox66, Sites.loxP, Sites.lox72),

            Excision(Sites.loxP, Sites.loxP, Sites.loxP),
            Excision(Sites.lox66, Sites.loxP, Sites.loxP),
            Excision(Sites.loxP, Sites.lox66, Sites.lox66),
            Excision(Sites.lox71, Sites.loxP, Sites.lox71),
            Excision(Sites.loxP, Sites.lox71, Sites.loxP),

            Inversion(Sites.loxP, rc(Sites.loxP), Sites.loxP, rc(Sites.loxP)),
            Inversion(Sites.lox66, rc(Sites.loxP), Sites.loxP, rc(Sites.lox66)),
            Inversion(Sites.lox71, rc(Sites.loxP), Sites.lox71, rc(Sites.loxP)),
        ]
