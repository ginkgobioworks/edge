from edge.ssr import lower_no_whitespace, rc
from edge.tests.test_ssr_view import SSRTester


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


class CreLoxTest(SSRTester):

    def test_integration_with_lox66_onto_lox71(self):
        self.with_genome("t" * 100 + lox71 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox66 + "g" * 10) \
            .modified_sequence_is("t" * 100 + lox72 + "g" * 10 + "a" * 10 + loxP + "c" * 100)

    def test_integration_with_lox66_onto_rc_lox71(self):
        self.with_genome("t" * 100 + rc(lox71) + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox66 + "g" * 10) \
            .modified_sequence_is(
                "t" * 100 + rc(loxP) + "t" * 10 + "c" * 10 + rc(lox72) + "c" * 100
            )

    def test_integration_with_lox71_onto_lox66(self):
        self.with_genome("t" * 100 + lox66 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox71 + "g" * 10) \
            .modified_sequence_is("t" * 100 + loxP + "g" * 10 + "a" * 10 + lox72 + "c" * 100)

    def test_integration_with_lox71_onto_rc_lox66(self):
        self.with_genome("t" * 100 + rc(lox66) + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox71 + "g" * 10) \
            .modified_sequence_is(
                "t" * 100 + rc(lox72) + "t" * 10 + "c" * 10 + rc(loxP) + "c" * 100
            )

    def test_integration_with_lox66_onto_loxP_is_not_allowed_because_immediate_loopout(self):
        self.with_genome("t" * 100 + loxP + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox66 + "g" * 10) \
            .is_not_allowed()

    def test_excision_with_loxP(self):
        self.with_genome("t" * 100 + loxP + "aaa" + loxP + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "c" * 100)

    def test_excision_with_lox66_loxP(self):
        self.with_genome("t" * 100 + lox66 + "aaa" + loxP + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "c" * 100)

    def test_excision_with_loxP_lox66(self):
        self.with_genome("t" * 100 + loxP + "aaa" + lox66 + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + lox66 + "c" * 100)

    def test_excision_with_lox71_loxP(self):
        self.with_genome("t" * 100 + lox71 + "aaa" + loxP + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + lox71 + "c" * 100)

    def test_excision_with_loxP_lox71(self):
        self.with_genome("t" * 100 + loxP + "aaa" + lox71 + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "c" * 100)

    def test_inversion_with_loxP(self):
        self.with_genome("t" * 100 + loxP + "atg" + rc(loxP) + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "cat" + rc(loxP) + "c" * 100)

    def test_inversion_with_loxP_lox66(self):
        self.with_genome("t" * 100 + loxP + "atg" + rc(lox66) + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + lox66 + "cat" + rc(loxP) + "c" * 100)

    def test_inversion_with_loxP_lox71(self):
        self.with_genome("t" * 100 + loxP + "atg" + rc(lox71) + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "cat" + rc(lox71) + "c" * 100)

    def test_inversion_with_lox66_loxP(self):
        self.with_genome("t" * 100 + lox66 + "atg" + rc(loxP) + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + loxP + "cat" + rc(lox66) + "c" * 100)

    def test_inversion_with_lox71_loxP(self):
        self.with_genome("t" * 100 + lox71 + "atg" + rc(loxP) + "c" * 100) \
            .when_trigger_crelox() \
            .modified_sequence_is("t" * 100 + lox71 + "cat" + rc(loxP) + "c" * 100)

    def test_RMCE_with_lox2272_loxP(self):
        self.with_genome("t" * 100 + lox2272 + "atg" + loxP + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox2272 + "g" * 10 + loxP, False) \
            .modified_sequence_is("t" * 100 + lox2272 + "g" * 10 + loxP + "c" * 100)

    def test_RMCE_with_loxP_lox2272(self):
        self.with_genome("t" * 100 + loxP + "atg" + lox2272 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + loxP + "g" * 10 + lox2272, False) \
            .modified_sequence_is("t" * 100 + loxP + "g" * 10 + lox2272 + "c" * 100)

    def test_RMCE_with_lox2272_loxP_and_reversed_on_donor(self):
        self.with_genome("t" * 100 + lox2272 + "atg" + loxP + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + rc(loxP) + "g" * 10 + rc(lox2272), False) \
            .modified_sequence_is("t" * 100 + lox2272 + "c" * 10 + loxP + "c" * 100)

    def test_RMCE_with_lox2272_loxP_and_loxP_lox2272_is_not_allowed(self):
        self.with_genome("t" * 100 + lox2272 + "atg" + loxP + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + loxP + "g" * 10 + lox2272, False) \
            .is_not_allowed()

    def test_RMCE_with_lox71_lox2272_and_lox66_lox2272(self):
        # this also tests that RMCE takes precedence over integration
        self.with_genome("t" * 100 + lox71 + "atg" + lox2272 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox66 + "g" * 10 + lox2272, False) \
            .modified_sequence_is("t" * 100 + lox72 + "g" * 10 + lox2272 + "c" * 100)

    def test_RMCE_with_lox2272_lox71_and_lox2272_lox66(self):
        # this also tests that RMCE takes precedence over integration
        self.with_genome("t" * 100 + lox2272 + "atg" + lox71 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + lox2272 + "g" * 10 + lox66, False) \
            .modified_sequence_is("t" * 100 + lox2272 + "g" * 10 + lox72 + "c" * 100)

    def test_RMCE_with_loxm2_66_lox71_and_loxm2_71_lox66(self):
        # this also tests that RMCE takes precedence over integration
        self.with_genome("t" * 100 + loxm2_66 + "atg" + lox71 + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + loxm2_71 + "g" * 10 + lox66) \
            .modified_sequence_is("t" * 100 + loxP + "g" * 10 + lox72 + "c" * 100)

    def test_RMCE_with_loxP_loxP_and_loxP_loxP_is_not_allowed_because_multiple_recombinations(self):
        self.with_genome("t" * 100 + loxP + "atg" + loxP + "c" * 100) \
            .when_trigger_crelox_with_donor("a" * 10 + loxP + "g" * 10 + loxP) \
            .is_not_allowed()
