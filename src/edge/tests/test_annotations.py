from django.test import TestCase

from edge.models import (
    Genome,
    Fragment,
    Operation,
)


class AnnotationsTest(TestCase):

    def setUp(self):
        self.genome = Genome.create('Test')
        self.root_sequence = 'agttcgaggctga'
        self.root = Fragment.create_with_sequence('Foo', self.root_sequence)

    def test_annotate_whole_chunk(self):
        self.root.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 1)
        self.assertEquals(self.root.annotations()[0].base_first, 1)
        self.assertEquals(self.root.annotations()[0].base_last, len(self.root_sequence))
        self.assertEquals(self.root.annotations()[0].feature.name, 'A1')
        self.assertEquals(self.root.annotations()[0].feature_base_first, 1)
        self.assertEquals(self.root.annotations()[0].feature_base_last, len(self.root_sequence))

    def test_annotate_part_of_chunk(self):
        self.root.annotate(1, 3, 'A1', 'gene', 1)
        self.root.annotate(7, 9, 'A2', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 2)
        self.assertEquals(self.root.annotations()[0].base_first, 1)
        self.assertEquals(self.root.annotations()[0].base_last, 3)
        self.assertEquals(self.root.annotations()[0].feature.name, 'A1')
        self.assertEquals(self.root.annotations()[0].feature_base_first, 1)
        self.assertEquals(self.root.annotations()[0].feature_base_last, 3)
        self.assertEquals(self.root.annotations()[1].base_first, 7)
        self.assertEquals(self.root.annotations()[1].base_last, 9)
        self.assertEquals(self.root.annotations()[1].feature.name, 'A2')
        self.assertEquals(self.root.annotations()[1].feature_base_first, 1)
        self.assertEquals(self.root.annotations()[1].feature_base_last, 3)

    def test_annotate_with_operation(self):
        op = Operation(genome=self.genome, type=Operation.RECOMBINATION[0])
        op.save()
        self.root.annotate(1, 3, 'A1', 'gene', 1, operation=op)
        self.assertEquals(len(self.root.annotations()), 1)
        self.assertEquals(self.root.annotations()[0].feature.operation, op)

    def test_annotate_multiple_chunks(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        # fragment now has three chunks
        u.annotate(2, 8, 'A1', 'gene', 1)
        self.assertEquals(len(u.annotations()), 1)
        self.assertEquals(u.annotations()[0].base_first, 2)
        self.assertEquals(u.annotations()[0].base_last, 8)
        self.assertEquals(u.annotations()[0].feature.name, 'A1')
        self.assertEquals(u.annotations()[0].feature_base_first, 1)
        self.assertEquals(u.annotations()[0].feature_base_last, 7)

    def test_annotation_adjusted_after_insert_outside_of_annotation(self):
        self.root.annotate(7, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 13)
        self.assertEquals(f.annotations()[0].base_last, 15)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_remove_outside_of_annotation(self):
        self.root.annotate(7, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.remove_bases(3, 4)
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 3)
        self.assertEquals(f.annotations()[0].base_last, 5)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_replace_outside_of_annotation(self):
        self.root.annotate(7, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.replace_bases(3, 4, 'cc')
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 5)
        self.assertEquals(f.annotations()[0].base_last, 7)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_annotation_adjusted_after_insert_within_annotation(self):
        self.root.annotate(2, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
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
        self.root.annotate(2, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.remove_bases(3, 4)
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
        self.root.annotate(2, 9, 'A1', 'gene', 1)
        f = self.root.update('Bar')
        f.replace_bases(3, 4, 'cc')
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
        f.annotate(9, 3, 'A1', 'gene', 1)
        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 1)
        self.assertEquals(f.annotations()[0].base_last, 3)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first,
                          len(self.root_sequence) - 9 + 1 + 1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence) - 9 + 1 + 3)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence) - 9 + 1 + 3)
        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 1)
        self.assertEquals(f.annotations()[1].feature_base_last, len(self.root_sequence) - 9 + 1)
        self.assertEquals(f.annotations()[1].feature.length, len(self.root_sequence) - 9 + 1 + 3)

    def test_annotate_circular_fragment_ending_at_bp_1(self):
        f = Fragment.create_with_sequence('Foo', self.root_sequence, circular=True)
        f.annotate(9, 1, 'A1', 'gene', 1)
        self.assertEquals(len(f.annotations()), 2)
        self.assertEquals(f.annotations()[0].base_first, 1)
        self.assertEquals(f.annotations()[0].base_last, 1)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first,
                          len(self.root_sequence) - 9 + 1 + 1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence) - 9 + 1 + 1)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence) - 9 + 1 + 1)
        self.assertEquals(f.annotations()[1].base_first, 9)
        self.assertEquals(f.annotations()[1].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[1].feature.name, 'A1')
        self.assertEquals(f.annotations()[1].feature_base_first, 1)
        self.assertEquals(f.annotations()[1].feature_base_last, len(self.root_sequence) - 9 + 1)
        self.assertEquals(f.annotations()[1].feature.length, len(self.root_sequence) - 9 + 1 + 1)

    def test_annotate_circular_fragment_ending_at_base_last(self):
        f = Fragment.create_with_sequence('Foo', self.root_sequence, circular=True)
        f.annotate(9, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 9)
        self.assertEquals(f.annotations()[0].base_last, len(self.root_sequence))
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, len(self.root_sequence) - 9 + 1)
        self.assertEquals(f.annotations()[0].feature.length, len(self.root_sequence) - 9 + 1)

    def test_annotations_list_ordered_by_bp(self):
        # annotate bps 5-9 before bps 1-6, but ordered annotations should still
        # return annotation of 1-6 before annotation of 5-9
        self.root.annotate(5, 9, 'A2', 'gene', 1)
        self.root.annotate(1, 6, 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 2)
        self.assertEquals(self.root.annotations()[0].base_first, 1)
        self.assertEquals(self.root.annotations()[0].base_last, 6)
        self.assertEquals(self.root.annotations()[0].feature.name, 'A1')
        self.assertEquals(self.root.annotations()[0].feature_base_first, 1)
        self.assertEquals(self.root.annotations()[0].feature_base_last, 6)
        self.assertEquals(self.root.annotations()[1].base_first, 5)
        self.assertEquals(self.root.annotations()[1].base_last, 9)
        self.assertEquals(self.root.annotations()[1].feature.name, 'A2')
        self.assertEquals(self.root.annotations()[1].feature_base_first, 1)
        self.assertEquals(self.root.annotations()[1].feature_base_last, 5)

    def test_getting_annotations_by_bp(self):
        self.root.annotate(1, 3, 'A1', 'gene', 1)
        self.root.annotate(7, 9, 'A2', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 2)
        self.assertEquals(len(self.root.annotations(bp_lo=1, bp_hi=3)), 1)
        self.assertEquals(self.root.annotations(bp_lo=1, bp_hi=3)[0].base_first, 1)
        self.assertEquals(self.root.annotations(bp_lo=1, bp_hi=3)[0].base_last, 3)
        self.assertEquals(self.root.annotations(bp_lo=1, bp_hi=3)[0].feature.name, 'A1')
        self.assertEquals(self.root.annotations(bp_lo=1, bp_hi=3)[0].feature_base_first, 1)
        self.assertEquals(self.root.annotations(bp_lo=1, bp_hi=3)[0].feature_base_last, 3)
        self.assertEquals(len(self.root.annotations(bp_lo=4, bp_hi=9)), 1)
        self.assertEquals(self.root.annotations(bp_lo=4, bp_hi=9)[0].base_first, 7)
        self.assertEquals(self.root.annotations(bp_lo=4, bp_hi=9)[0].base_last, 9)
        self.assertEquals(self.root.annotations(bp_lo=4, bp_hi=9)[0].feature.name, 'A2')
        self.assertEquals(self.root.annotations(bp_lo=4, bp_hi=9)[0].feature_base_first, 1)
        self.assertEquals(self.root.annotations(bp_lo=4, bp_hi=9)[0].feature_base_last, 3)

    def test_getting_annotations_by_bp_not_at_chunk_break_points(self):
        self.root.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations(bp_lo=1, bp_hi=len(self.root_sequence))), 1)
        self.assertEquals(len(self.root.annotations(bp_lo=2, bp_hi=8)), 1)

    def test_merge_splitted_annotations(self):
        self.root.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 1)

        f = self.root.update('Bar')
        f._find_and_split_before(4)

        n = len(self.root_sequence)
        self.assertEquals(len(f.annotations(bp_lo=1, bp_hi=n)), 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].base_last, n)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature.name, 'A1')
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_first, 1)
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[0].feature_base_last, n)

    def test_does_not_merge_splitted_annotations_if_bp_removed_from_annotation(self):
        self.root.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 1)

        f = self.root.update('Bar')
        f.remove_bases(4, 3)

        n = len(self.root_sequence) - 3
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
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_last, n + 3)

    def test_does_not_merge_splitted_annotations_if_bp_inserted_in_annotation(self):
        self.root.annotate(1, len(self.root_sequence), 'A1', 'gene', 1)
        self.assertEquals(len(self.root.annotations()), 1)

        f = self.root.update('Bar')
        f.insert_bases(4, 'ccc')

        n = len(self.root_sequence) + 3
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
        self.assertEquals(f.annotations(bp_lo=1, bp_hi=n)[1].feature_base_last, n - 3)

    def test_child_inherits_annotation_from_parent_when_parent_annotates_preserved_region(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')

        self.assertEquals(len(f.annotations()), 0)
        self.root.annotate(10, 12, 'A1', 'gene', 1)

        f = Fragment.objects.get(pk=f.pk).indexed_fragment()
        self.assertEquals(len(f.annotations()), 1)
        self.assertEquals(f.annotations()[0].base_first, 16)
        self.assertEquals(f.annotations()[0].base_last, 18)
        self.assertEquals(f.annotations()[0].feature.name, 'A1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)

    def test_child_inherits_annotation_from_parent_when_parent_annotates_updated_region(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')

        self.assertEquals(len(f.annotations()), 0)
        self.root.annotate(2, 9, 'A1', 'gene', 1)

        f = Fragment.objects.get(pk=f.pk).indexed_fragment()
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
        new_f.annotate(2, 4, 'X1', 'gene', 1)

        self.assertEquals(len(new_f.annotations()), 1)
        self.assertEquals(len(self.root.annotations()), 0)

        f = self.root.update('Bar')
        f.insert_fragment(3, new_f)

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

        f = self.root.update('Bar')
        f.insert_fragment(3, new_f)

        self.assertEquals(len(new_f.annotations()), 0)
        self.assertEquals(len(self.root.annotations()), 0)
        self.assertEquals(len(f.annotations()), 0)

        new_f.annotate(2, 4, 'X1', 'gene', 1)

        self.root = Fragment.objects.get(pk=self.root.pk).indexed_fragment()
        f = Fragment.objects.get(pk=f.pk).indexed_fragment()
        self.assertEquals(len(new_f.annotations()), 1)
        self.assertEquals(len(self.root.annotations()), 0)
        self.assertEquals(len(f.annotations()), 1)

        self.assertEquals(f.annotations()[0].base_first, 4)
        self.assertEquals(f.annotations()[0].base_last, 6)
        self.assertEquals(f.annotations()[0].feature.name, 'X1')
        self.assertEquals(f.annotations()[0].feature_base_first, 1)
        self.assertEquals(f.annotations()[0].feature_base_last, 3)
