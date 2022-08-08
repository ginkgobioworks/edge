from edge.ssr import lower_no_whitespace, rc
from edge.tests.test_ssr_view import SSRTester


Frt = lower_no_whitespace("GAAGTTCCTATTC tctagaaa GTATAGGAACTTC")
Frt3 = lower_no_whitespace("GAAGTTCCTATTC ttcaaata GTATAGGAACTTC")


class FlpTest(SSRTester):

    def test_excise_with_frt(self):
        self.with_genome("t" * 100 + Frt + "aaa" + Frt + "c" * 100) \
            .when_trigger_flp() \
            .modified_sequence_is("t" * 100 + Frt + "c" * 100)

        self.with_genome("t" * 100 + Frt3 + "aaa" + Frt3 + "c" * 100) \
            .when_trigger_flp() \
            .modified_sequence_is("t" * 100 + Frt3 + "c" * 100)

    def test_inversion_with_frt(self):
        self.with_genome("t" * 100 + Frt + "atg" + rc(Frt) + "c" * 100) \
            .when_trigger_flp() \
            .modified_sequence_is("t" * 100 + Frt + "cat" + rc(Frt) + "c" * 100)

        self.with_genome("t" * 100 + Frt3 + "atg" + rc(Frt3) + "c" * 100) \
            .when_trigger_flp() \
            .modified_sequence_is("t" * 100 + Frt3 + "cat" + rc(Frt3) + "c" * 100)

    def test_RMCE_with_frt_frt3(self):
        self.with_genome("t" * 100 + Frt + "atg" + Frt3 + "c" * 100) \
            .when_trigger_flp_with_donor("a" * 10 + Frt + "g" * 10 + Frt3, False) \
            .modified_sequence_is("t" * 100 + Frt + "g" * 10 + Frt3 + "c" * 100)

    def test_RMCE_with_frt3_frt(self):
        self.with_genome("t" * 100 + Frt3 + "atg" + Frt + "c" * 100) \
            .when_trigger_flp_with_donor("a" * 10 + Frt3 + "g" * 10 + Frt, False) \
            .modified_sequence_is("t" * 100 + Frt3 + "g" * 10 + Frt + "c" * 100)
