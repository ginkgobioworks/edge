from django.test import TestCase
from edge.models import *


class AnnotationsTest(TestCase):

    def setUp(self):
        self.root_sequence = 'agttcgaggctga'
        self.root = Fragment.create_with_sequence('Foo', self.root_sequence)

    def test_annotate_whole_chunk(self):
        a = self.root.annotate()
        a.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(a.annotations()), 1)
        self.assertEquals(a.annotations()[0].base_first, 1)
        self.assertEquals(a.annotations()[0].base_last, len(self.root_sequence))
        self.assertEquals(a.annotations()[0].feature.name, 'A1')
        self.assertEquals(a.annotations()[0].feature_base_first, 1)
        self.assertEquals(a.annotations()[0].feature_base_last, len(self.root_sequence))

    def test_annotate_part_of_chunk(self):
        a = self.root.annotate()
        a.annotate(1, 3, 'A1', 'gene', 1)
        a.annotate(7, 9, 'A2', 'gene', 1)
        self.assertEquals(len(a.annotations()), 2)
        self.assertEquals(a.annotations()[0].base_first, 1)
        self.assertEquals(a.annotations()[0].base_last, 3)
        self.assertEquals(a.annotations()[0].feature.name, 'A1')
        self.assertEquals(a.annotations()[0].feature_base_first, 1)
        self.assertEquals(a.annotations()[0].feature_base_last, 3)
        self.assertEquals(a.annotations()[1].base_first, 7)
        self.assertEquals(a.annotations()[1].base_last, 9)
        self.assertEquals(a.annotations()[1].feature.name, 'A2')
        self.assertEquals(a.annotations()[1].feature_base_first, 1)
        self.assertEquals(a.annotations()[1].feature_base_last, 3)

    def test_annotate_multiple_chunks(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u
        # f now has three chunks
        a = f.annotate()
        a.annotate(2, 8, 'A1', 'gene', 1)
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 2)
        self.assertEquals(f.annotations()[0].base_last, 8)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 7)

    def test_annotation_adjusted_after_insert_outside_of_annotation(self):
        a = self.root.annotate()
        a.annotate(7, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 13)
        self.assertEquals(f.annotations()[0].base_last, 15)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_remove_outside_of_annotation(self):
        a = self.root.annotate()
        a.annotate(7, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.remove_bases(3, 4)
        f = u
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 3)
        self.assertEquals(f.annotations()[0].base_last, 5)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_replace_outside_of_annotation(self):
        a = self.root.annotate()
        a.annotate(7, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.replace_bases(3, 4, 'cc')
        f = u
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 5)
        self.assertEquals(f.annotations()[0].base_last, 7)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_insert_within_annotation(self):
        a = self.root.annotate()
        a.annotate(2, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u
        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 2)
        self.assertEquals(f.annotations()[0].base_last, 2)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 1)
        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, 15)
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 2)
        self.assertEquals(f.annotations()[1].feature_base_last, 8)

    def test_annotation_adjusted_after_remove_within_annotation(self):
        a = self.root.annotate()
        a.annotate(2, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.remove_bases(3, 4)
        f = u
        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 2)
        self.assertEquals(f.annotations()[0].base_last, 2)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 1)
        self.assertEquals(f.annotations()[1].base_first, 3)
        self.assertEquals(f.annotations()[1].base_last, 5)
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 6)
        self.assertEquals(f.annotations()[1].feature_base_last, 8)

    def test_annotation_adjusted_after_replace_within_annotation(self):
        a = self.root.annotate()
        a.annotate(2, 9, 'A1', 'gene', 1)
        u = self.root.update('Bar')
        u.replace_bases(3, 4, 'cc')
        f = u
        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 2)
        self.assertEquals(f.annotations()[0].base_last, 2)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 1)
        self.assertEquals(f.annotations()[1].base_first, 5)
        self.assertEquals(f.annotations()[1].base_last, 7)
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 6)
        self.assertEquals(f.annotations()[1].feature_base_last, 8)

    def test_annotate_circular_fragment(self):
        f = Fragment.create_with_sequence('Foo', self.root_sequence, circular=True)
        a = f.annotate()
        a.annotate(9, 3, 'A1', 'gene', 1)

        self.assertEquals(len(f.annotations()), 2)

        self.assertEquals(f.annotations()[0].base_first, 1)
        self.assertEquals(f.annotations()[0].base_last, 3)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, len(self.root_sequence)-9+1+1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence)-9+1+3)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence)-9+1+3)

        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 1)
        self.assertEquals(f.annotations()[1].feature_base_last, len(self.root_sequence)-9+1)
        self.assertEquals(f.annotations()[1].feature.length, len(self.root_sequence)-9+1+3)

    def test_annotate_circular_fragment_ending_at_bp_1(self):
        f = Fragment.create_with_sequence('Foo', self.root_sequence, circular=True)
        a = f.annotate()
        a.annotate(9, 1, 'A1', 'gene', 1)

        self.assertEquals(len(f.annotations()), 2)

        self.assertEquals(f.annotations()[0].base_first, 1)
        self.assertEquals(f.annotations()[0].base_last, 1)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, len(self.root_sequence)-9+1+1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence)-9+1+1)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence)-9+1+1)

        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 1)
        self.assertEquals(f.annotations()[1].feature_base_last, len(self.root_sequence)-9+1)
        self.assertEquals(f.annotations()[1].feature.length, len(self.root_sequence)-9+1+1)

    def test_annotate_circular_fragment_ending_at_base_last(self):
        f = Fragment.create_with_sequence('Foo', self.root_sequence, circular=True)
        a = f.annotate()
        a.annotate(9, len(self.root_sequence), 'A1', 'gene', 1)

        self.assertEquals(len(f.annotations()), 1)

        self.assertEquals(f.annotations()[0].base_first, 9)
        self.assertEquals(f.annotations()[0].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence)-9+1)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence)-9+1)

    def test_getting_annotations_by_bp(self):
        a = self.root.annotate()
        a.annotate(1, 3, 'A1', 'gene', 1)
        a.annotate(7, 9, 'A2', 'gene', 1)

        self.assertEquals(len(a.annotations()), 2)

        self.assertEquals(len(a.annotations(bp_lo=1, bp_hi=3)), 1)
        self.assertEquals(a.annotations(bp_lo=1, bp_hi=3)[0].base_first, 1)
        self.assertEquals(a.annotations(bp_lo=1, bp_hi=3)[0].base_last, 3)
        self.assertEquals(a.annotations(bp_lo=1, bp_hi=3)[0].feature.name, 'A1')
        self.assertEquals(a.annotations(bp_lo=1, bp_hi=3)[0].feature_base_first, 1)
        self.assertEquals(a.annotations(bp_lo=1, bp_hi=3)[0].feature_base_last, 3)

        self.assertEquals(len(a.annotations(bp_lo=4, bp_hi=9)), 1)
        self.assertEquals(a.annotations(bp_lo=4, bp_hi=9)[0].base_first, 7)
        self.assertEquals(a.annotations(bp_lo=4, bp_hi=9)[0].base_last, 9)
        self.assertEquals(a.annotations(bp_lo=4, bp_hi=9)[0].feature.name, 'A2')
        self.assertEquals(a.annotations(bp_lo=4, bp_hi=9)[0].feature_base_first, 1)
        self.assertEquals(a.annotations(bp_lo=4, bp_hi=9)[0].feature_base_last, 3)

    def test_getting_annotations_by_bp_not_at_chunk_break_points(self):
        a = self.root.annotate()
        a.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(a.annotations(bp_lo=1, bp_hi=len(self.root_sequence))), 1)
        self.assertEquals(len(a.annotations(bp_lo=2, bp_hi=8)), 1)

    def test_merge_splitted_annotations(self):
        a = self.root.annotate()
        a.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(a.annotations()), 1)

        u = self.root.update('Bar')
        u._find_and_split_before(4)
        f = u
        n = len(self.root_sequence)
        self.assertEquals(len(f.annotations(bp_lo=1, bp_hi=n)), 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_last, n)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_last, n)

    def test_does_not_merge_splitted_annotations_if_bp_removed_from_annotation(self):
        a = self.root.annotate()
        a.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(a.annotations()), 1)

        u = self.root.update('Bar')
        u.remove_bases(4, 3)
        f = u
        n = len(self.root_sequence)-3
        self.assertEquals(len(f.annotations(bp_lo=1, bp_hi=n)), 2)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_last, 3)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_last, 3)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].base_first, 4)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].base_last, n)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_first, 7)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_last, n+3)

    def test_does_not_merge_splitted_annotations_if_bp_inserted_in_annotation(self):
        a = self.root.annotate()
        a.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(a.annotations()), 1)

        u = self.root.update('Bar')
        u.insert_bases(4, 'ccc')
        f = u
        n = len(self.root_sequence)+3
        self.assertEquals(len(f.annotations(bp_lo=1, bp_hi=n)), 2)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_last, 3)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_last, 3)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].base_first, 7)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].base_last, n)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_first, 4)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_last, n-3)

    def test_child_inherits_annotation_from_parent_when_parent_annotates_preserved_region(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')

        f = u
        self.assertEquals(len(f.annotations()), 0)
        a = self.root.annotate()
        a.annotate(10, 12, 'A1', 'gene', 1)

        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 16)
        self.assertEquals(f.annotations()[0].base_last, 18)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_child_inherits_annotation_from_parent_when_parent_annotates_updated_region(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u

        self.assertEquals(len(f.annotations()), 0)
        a = self.root.annotate()
        a.annotate(2, 9, 'A1', 'gene', 1)

        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 2)
        self.assertEquals(f.annotations()[0].base_last, 2)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 1)
        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, 15)
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 2)
        self.assertEquals(f.annotations()[1].feature_base_last, 8)

    def test_inherits_annotations_from_inserted_fragment(self):
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        a = new_f.annotate()
        a.annotate(2, 4, 'X1', 'gene', 1)

        self.assertEquals(len(new_f.annotations()), 1)
        self.assertEquals(len(self.root.annotations()), 0)

        u = self.root.update('Bar')
        u.insert_fragment(3, new_f)
        f = u

        self.assertEquals(len(new_f.annotations()), 1)
        self.assertEquals(len(self.root.annotations()), 0)
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 4)
        self.assertEquals(f.annotations()[0].base_last, 6)
        self.assertEquals(f.annotations()[0].feature.name, 'X1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_inherits_new_annotations_on_inserted_fragment(self):
        new_f = Fragment.create_with_sequence('Test', 'gataca')

        self.assertEquals(len(new_f.annotations()), 0)
        self.assertEquals(len(self.root.annotations()), 0)

        u = self.root.update('Bar')
        u.insert_fragment(3, new_f)
        f = u

        self.assertEquals(len(new_f.annotations()), 0)
        self.assertEquals(len(self.root.annotations()), 0)
        self.assertEquals(len(f.annotations()), 0)

        a = new_f.annotate()
        a.annotate(2, 4, 'X1', 'gene', 1)

        self.assertEquals(len(new_f.annotations()), 1)
        self.assertEquals(len(self.root.annotations()), 0)
        self.assertEquals(len(f.annotations()), 1)

        self.assertEquals(f.annotations()[0].base_first, 4)
        self.assertEquals(f.annotations()[0].base_last, 6)
        self.assertEquals(f.annotations()[0].feature.name, 'X1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)
