from django.test import TestCase
from edge.models import *


class UserDefinedFragmentsTests(TestCase):

    def test_user_defined_fragments_does_not_include_genomic_fragment(self):
        genome = Genome.create('Foo')
        s = 'atggcatattcgcagct'
        f1 = genome.add_fragment('chrI', s)

        f2 = Fragment.create_with_sequence('Bar', 'aacctaaaattataa')
        self.assertEquals(len(Fragment.user_defined_fragments()), 1)
        self.assertEquals(Fragment.user_defined_fragments()[0].name, 'Bar')
        self.assertEquals(Fragment.user_defined_fragments()[0].id, f2.id)

    def test_user_defined_fragments_does_not_include_inactive_fragments(self):
        f1 = Fragment.create_with_sequence('Bar', 'aacctaaaattataa')
        self.assertEquals(len(Fragment.user_defined_fragments()), 1)
        self.assertEquals(Fragment.user_defined_fragments()[0].id, f1.id)
        f1.active = False
        f1.save()
        self.assertEquals(len(Fragment.user_defined_fragments()), 0)
