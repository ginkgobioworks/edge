import tempfile

from django.test import TestCase
from edge.models import Fragment, Genome


class UserDefinedFragmentsTests(TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_user_defined_fragments_does_not_include_genomic_fragment(self):
        genome = Genome.create("Foo")
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s, dirn=self.tmpdir.name)

        f = Fragment.create_with_sequence(
            "Bar", "aacctaaaattataa", dirn=self.tmpdir.name
        )
        self.assertEquals(len(Fragment.user_defined_fragments()), 1)
        self.assertEquals(Fragment.user_defined_fragments()[0].name, "Bar")
        self.assertEquals(Fragment.user_defined_fragments()[0].id, f.id)

    def test_user_defined_fragments_does_not_include_inactive_fragments(self):
        f1 = Fragment.create_with_sequence(
            "Bar", "aacctaaaattataa", dirn=self.tmpdir.name
        )
        self.assertEquals(len(Fragment.user_defined_fragments()), 1)
        self.assertEquals(Fragment.user_defined_fragments()[0].id, f1.id)
        f1.active = False
        f1.save()
        self.assertEquals(len(Fragment.user_defined_fragments()), 0)
