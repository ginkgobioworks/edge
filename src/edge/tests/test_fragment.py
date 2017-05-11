from django.test import TestCase

from edge.models import (
    Chunk,
    Edge,
    Fragment,
)


class FragmentCreateTests(TestCase):

    def test_can_create_fragment_with_no_chunk_size(self):
        f = Fragment.create_with_sequence('Foo', 'gataccggtactag', initial_chunk_size=None)
        self.assertEquals(f.sequence, 'gataccggtactag')

    def test_can_create_fragment_with_different_chunk_sizes(self):
        s = 'gataccggtactag'
        f = Fragment.create_with_sequence('Foo', s, initial_chunk_size=len(s))
        self.assertEquals(f.sequence, s)
        f = Fragment.create_with_sequence('Foo', s, initial_chunk_size=0)
        self.assertEquals(f.sequence, s)
        f = Fragment.create_with_sequence('Foo', s, initial_chunk_size=1)
        self.assertEquals(f.sequence, s)
        f = Fragment.create_with_sequence('Foo', s, initial_chunk_size=3)
        self.assertEquals(f.sequence, s)
        f = Fragment.create_with_sequence('Foo', s, initial_chunk_size=len(s) * 1000)
        self.assertEquals(f.sequence, s)


class FragmentTests(TestCase):

    def setUp(self):
        self.root_sequence = 'agttcgaggctga'
        self.root = Fragment.create_with_sequence('Foo', self.root_sequence)

    def test_root_fragment(self):
        self.assertEquals(self.root.sequence, self.root_sequence)
        self.assertEquals(self.root.name, 'Foo')

    def test_split_in_middle(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(4)
        self.assertEquals([c.sequence for c in u.chunks()], ['agt', 'tcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_can_reconstruct_index_after_split(self):
        u = self.root.update('Bar')
        u._find_and_split_before(4)
        self.assertEquals(u.sequence, self.root_sequence)
        self.assertEquals(u.length, len(self.root_sequence))
        u.fragment_chunk_location_set.all().delete()
        u.index_fragment_chunk_locations()
        self.assertEquals(u.sequence, self.root_sequence)
        self.assertEquals(u.length, len(self.root_sequence))

    def test_split_at_start(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(1)
        self.assertEquals([c.sequence for c in u.chunks()], ['agttcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_split_at_end(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(len(self.root_sequence))
        self.assertEquals([c.sequence for c in u.chunks()], ['agttcgaggctg', 'a'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_split_past_end(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(len(self.root_sequence) + 1)
        self.assertEquals([c.sequence for c in u.chunks()], ['agttcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_double_split_at_same_position(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(3)
        self.assertEquals([c.sequence for c in u.chunks()], ['ag', 'ttcgaggctga'])
        prev, cur = u._find_and_split_before(3)
        self.assertEquals([c.sequence for c in u.chunks()], ['ag', 'ttcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_split_chunk_with_size_of_one(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(4)
        self.assertEquals([c.sequence for c in u.chunks()], ['agt', 'tcgaggctga'])
        prev, cur = u._find_and_split_before(3)
        self.assertEquals([c.sequence for c in u.chunks()], ['ag', 't', 'tcgaggctga'])
        prev, cur = u._find_and_split_before(4)
        self.assertEquals([c.sequence for c in u.chunks()], ['ag', 't', 'tcgaggctga'])
        prev, cur = u._find_and_split_before(3)
        self.assertEquals([c.sequence for c in u.chunks()], ['ag', 't', 'tcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

    def test_split_chunk_invalidates_location_indices(self):
        u1 = self.root.update('Bar')
        u1._find_and_split_before(4)
        u2 = self.root.update('Baz')

        self.assertEquals(self.root.has_location_index, True)
        self.assertEquals(u1.has_location_index, True)
        self.assertEquals(u2.has_location_index, True)

        # split, keeps parent (self.root) and current fragment (u2) index, but
        # invalidate indices of fragments using the splitted chunk
        u2._find_and_split_before(3)

        self.root = Fragment.objects.get(pk=self.root.pk)
        u1 = Fragment.objects.get(pk=u1.pk)
        u2 = Fragment.objects.get(pk=u2.pk)

        self.assertEquals(self.root.has_location_index, True)
        self.assertEquals(u1.has_location_index, False)
        self.assertEquals(u2.has_location_index, True)

        # can re-index u1
        u1 = u1.indexed_fragment()
        self.assertEquals(u1.sequence, self.root_sequence)
        self.assertEquals(u1.length, len(self.root_sequence))

    def test_split_does_not_invalidate_location_indices_if_not_splitting_a_chunk(self):
        u1 = self.root.update('Bar')
        u1._find_and_split_before(4)
        u2 = self.root.update('Baz')

        self.assertEquals(self.root.has_location_index, True)
        self.assertEquals(u1.has_location_index, True)
        self.assertEquals(u2.has_location_index, True)

        # not splitting a chunk, so location index of u1 not removed
        u2._find_and_split_before(4)

        self.assertEquals(self.root.has_location_index, True)
        self.assertEquals(u1.has_location_index, True)
        self.assertEquals(u2.has_location_index, True)

    def test_get_sequence_by_region(self):
        for i in range(0, len(self.root_sequence)):
            for j in range(i, len(self.root_sequence)):
                self.assertEquals(self.root.get_sequence(i + 1, j + 1), self.root_sequence[i:j + 1])

    def test_get_sequence_by_region_after_split(self):
        f = self.root.update('Bar')
        prev, cur = f._find_and_split_before(6)
        for i in range(0, len(self.root_sequence)):
            for j in range(i, len(self.root_sequence)):
                self.assertEquals(f.get_sequence(i + 1, j + 1), self.root_sequence[i:j + 1])

    def test_get_circular_sequence(self):
        s = 'agttcgaggctga'
        f = Fragment.create_with_sequence('Foo', s, circular=True)
        self.assertEquals(f.get_sequence(len(s) - 3 + 1, 4), s[-3:] + s[:4])

    def test_insert_sequence_in_middle(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2] + 'gataca' + self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_start(self):
        f = self.root.update('Bar')
        f.insert_bases(1, 'gataca')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'gataca' + self.root_sequence)
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_end(self):
        f = self.root.update('Bar')
        f.insert_bases(None, 'gataca')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence + 'gataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence + 'gatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence + 'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_does_not_allow_insert_same_fragment_twice_which_creates_loops(self):
        u = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        u.insert_fragment(3, new_f)
        self.assertRaises(Exception, u.insert_fragment, 3, new_f)

    def test_insert_new_fragment_in_middle(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.insert_fragment(3, new_f)
        self.assertEquals(
            f.sequence, self.root_sequence[0: 2] + 'gcccataca' + self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_start(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.insert_fragment(1, new_f)
        self.assertEquals(f.sequence, 'gcccataca' + self.root_sequence)
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_end(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.insert_fragment(None, new_f)
        self.assertEquals(f.sequence, self.root_sequence + 'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        u.insert_fragment(None, new_f)
        f1 = u
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence + 'gcccatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence + 'gcccataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_existing_fragment_in_middle(self):
        existing_f = Fragment.create_with_sequence('Test', 'gataca')
        f = self.root.update('Bar')
        f.insert_fragment(3, existing_f)
        self.assertEquals(f.sequence, self.root_sequence[0:2] + 'gataca' + self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_fragment_inherits_annotations_from_new_fragment(self):
        f = self.root.update('Bar')
        self.assertEquals([a for a in f.annotations() if a.name == 'Uma'], [])
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f.annotate(2, 4, 'Uma', 'feature', 1)
        f.insert_fragment(2, new_f)
        self.assertEquals([a for a in f.annotations() if a.feature.name == 'Uma'][0].base_first, 3)
        self.assertEquals([a for a in f.annotations() if a.feature.name == 'Uma'][0].base_last, 5)

    def test_remove_sequence_in_middle(self):
        f = self.root.update('Bar')
        f.remove_bases(3, 4)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2] + self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_start(self):
        f = self.root.update('Bar')
        f.remove_bases(1, 4)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[4:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_end(self):
        f = self.root.update('Bar')
        # removing before 4th to last bp, remove 4 bases
        f.remove_bases(len(self.root_sequence) - 3, 4)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_past_end(self):
        f = self.root.update('Bar')
        # removing before 2nd to last bp, remove 4 bases
        f.remove_bases(len(self.root_sequence) - 1, 4)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-2])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        # removing before 4th to last bp, remove 4 bases
        u.remove_bases(len(self.root_sequence) - 3, 4)
        f1 = u
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence[:-4] + 'aaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-4])
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_in_middle(self):
        f = self.root.update('Bar')
        f.replace_bases(3, 6, 'cccc')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2] + 'cccc' + self.root_sequence[8:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_start(self):
        f = self.root.update('Bar')
        f.replace_bases(1, 6, 'cccc')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'cccc' + self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_end(self):
        f = self.root.update('Bar')
        # removing before 6th to last bp, remove 6 bases
        f.replace_bases(len(self.root_sequence) - 5, 6, 'cccc')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-6] + 'cccc')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_past_end(self):
        f = self.root.update('Bar')
        # removing before 4th to last bp, remove 6 bases
        f.replace_bases(len(self.root_sequence) - 3, 6, 'cccc')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4] + 'cccc')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        # removing before 6th to last bp, remove 6 bases
        u.replace_bases(len(self.root_sequence) - 5, 6, 'cccc')
        f1 = u
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence[:-6] + 'ccccaaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-6] + 'cccc')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_in_middle(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.replace_with_fragment(3, 6, new_f)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(
            f.sequence, self.root_sequence[0: 2] + 'gcccataca' + self.root_sequence[8:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_start(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.replace_with_fragment(1, 6, new_f)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'gcccataca' + self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_end(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.replace_with_fragment(len(self.root_sequence) - 5, 6, new_f)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-6] + 'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_past_end(self):
        f = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        f.replace_with_fragment(len(self.root_sequence) - 3, 6, new_f)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4] + 'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        new_f = Fragment.create_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        u.replace_with_fragment(len(self.root_sequence) - 5, 6, new_f)
        f1 = u
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence[:-6] + 'gcccatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-6] + 'gcccataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_multiple_operations(self):
        f = self.root.update('Bar')
        # removing before 4th to last bp, remove 6 bases
        f.replace_bases(len(self.root_sequence) - 3, 6, 'cccc')
        f.insert_bases(None, 'ggg')
        f.remove_bases(1, 3)
        self.assertEquals(f.name, 'Bar')
        s = self.root_sequence
        self.assertEquals(f.sequence, s[3:-4] + 'cccc' + 'ggg')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_on_chunk_size_of_one(self):
        f = self.root.update('Bar')
        # this op creates a chunk with size of 1, then G, then the rest of the sequence
        f.insert_bases(2, 'G')
        # this op has to insert on the chunk with size of 1
        f.insert_bases(1, 'C')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'C' + self.root_sequence[0:1] + 'G' + self.root_sequence[1:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_double_insert(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'G')
        f.insert_bases(3, 'C')
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2] + 'CG' + self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_double_remove(self):
        f = self.root.update('Bar')
        f.remove_bases(3, 4)
        f.remove_bases(3, 4)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2] + self.root_sequence[10:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_at_end_for_two_fragments(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u
        u = f1.update('Baz')
        u.insert_bases(None, 'atatat')
        f2 = u
        self.assertEquals(f2.sequence, self.root_sequence + 'gataca' + 'atatat')
        self.assertEquals(f1.sequence, self.root_sequence + 'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_then_insert(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u
        u = f1.update('Baz')
        u.remove_bases(len(self.root_sequence) + 1, 6)
        f2 = u
        u = f2.update('Far')
        u.insert_bases(None, 'atatat')
        f3 = u
        self.assertEquals(f3.sequence, self.root_sequence + 'atatat')
        self.assertEquals(f2.sequence, self.root_sequence)
        self.assertEquals(f1.sequence, self.root_sequence + 'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_add_edges_replaces_existing_edge_with_same_fragment_id(self):
        u = self.root.update('Bar')
        prev, cur = u._find_and_split_before(4)
        self.assertEquals(prev.out_edges.count(), 1)
        self.assertEquals(prev.out_edges.all()[0].to_chunk_id, cur.id)
        self.assertEquals(prev.out_edges.all()[0].fragment_id, self.root.id)

        # add new edge, does not overwrite old edge from parent fragment
        fake_c = Chunk(id=100, initial_fragment=self.root)
        fake_c.save()
        u._add_edges(prev, Edge(from_chunk=prev, fragment_id=u.id, to_chunk=fake_c))
        self.assertEquals(prev.out_edges.count(), 2)
        self.assertEquals(prev.out_edges.all()[0].to_chunk_id, cur.id)
        self.assertEquals(prev.out_edges.all()[0].fragment_id, self.root.id)
        self.assertEquals(prev.out_edges.all()[1].to_chunk_id, 100)
        self.assertEquals(prev.out_edges.all()[1].fragment_id, u.id)

        # add another edge with fragment id 101, does overwrite edge added last
        # time, with same fragment id
        fake_c = Chunk(id=102, initial_fragment=self.root)
        fake_c.save()
        u._add_edges(prev, Edge(from_chunk=prev, fragment_id=u.id, to_chunk=fake_c))
        self.assertEquals(prev.out_edges.count(), 2)
        self.assertEquals(prev.out_edges.all()[0].to_chunk_id, cur.id)
        self.assertEquals(prev.out_edges.all()[0].fragment_id, self.root.id)
        self.assertEquals(prev.out_edges.all()[1].to_chunk_id, 102)
        self.assertEquals(prev.out_edges.all()[1].fragment_id, u.id)

    def test_find_chunk_by_edge_and_indexing_fragment_chunk_locations(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        f.annotate(3, 5, 'A1', 'gene', 1)

        # there should now be 4 chunks for f
        self.assertEquals([c.sequence for c in f.chunks()],
                          [self.root_sequence[0:2], 'gat', 'aca', self.root_sequence[2:]])
        self.assertEquals(f.has_location_index, True)

        # now remove all fragment chunk indices
        f.fragment_chunk_location_set.all().delete()
        self.assertEquals(f.has_location_index, False)

        # can still get chunks by walking edges
        self.assertEquals([c.sequence for c in f.chunks_by_walking()],
                          [self.root_sequence[0:2], 'gat', 'aca', self.root_sequence[2:]])

        # re-index
        f.index_fragment_chunk_locations()
        self.assertEquals(f.has_location_index, True)
        self.assertEquals([c.sequence for c in f.chunks()],
                          [self.root_sequence[0:2], 'gat', 'aca', self.root_sequence[2:]])

    def test_indexing_sets_est_length(self):
        self.root.est_length = None
        self.root.save()
        self.assertEquals(self.root.est_length, None)
        self.root.index_fragment_chunk_locations()
        self.assertEquals(self.root.est_length, len(self.root_sequence))

    def test_find_chunk(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')

        # there should now be 3 chunks for f
        self.assertEquals([c.sequence for c in f.chunks()],
                          [self.root_sequence[0:2], 'gataca', self.root_sequence[2:]])
        chunk_ids = [c.id for c in f.chunks()]

        # get first bp
        prev_chunk, chunk, next_chunk, bases_visited = f._find_chunk_prev_next(1)
        self.assertEquals(prev_chunk, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk.id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get last bp in first chunk
        prev_chunk, chunk, next_chunk, bases_visited = f._find_chunk_prev_next(2)
        self.assertEquals(prev_chunk, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk.id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get first bp in second chunk
        prev_chunk, chunk, next_chunk, bases_visited = f._find_chunk_prev_next(3)
        self.assertEquals(prev_chunk.id, chunk_ids[0])
        self.assertEquals(chunk.id, chunk_ids[1])
        self.assertEquals(next_chunk.id, chunk_ids[2])
        self.assertEquals(bases_visited, 2 + 6)

        # get last bp in last chunk
        prev_chunk, chunk, next_chunk, bases_visited = \
            f._find_chunk_prev_next(6 + len(self.root_sequence))
        self.assertEquals(prev_chunk.id, chunk_ids[1])
        self.assertEquals(chunk.id, chunk_ids[2])
        self.assertEquals(next_chunk, None)
        self.assertEquals(bases_visited, 6 + len(self.root_sequence))

        # get bp past last chunk
        prev_chunk, chunk, next_chunk, bases_visited = f._find_chunk_prev_next(100)
        self.assertEquals(prev_chunk.id, chunk_ids[2])
        self.assertEquals(chunk, None)
        self.assertEquals(next_chunk, None)
        self.assertEquals(bases_visited, 6 + len(self.root_sequence))


class FragmentChunkTest(TestCase):

    def setUp(self):
        self.root_sequence = 'agttcgaggctga'
        self.root = Fragment.create_with_sequence('Foo', self.root_sequence)

    def test_next_chunk(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(f.fragment_chunk(chunks[0]).next_chunk.id, chunks[1].id)
        self.assertEquals(f.fragment_chunk(chunks[1]).next_chunk.id, chunks[2].id)
        self.assertEquals(f.fragment_chunk(chunks[2]).next_chunk, None)

    def test_location(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(f.fragment_chunk(chunks[0]).location, (1, 2))
        self.assertEquals(f.fragment_chunk(chunks[1]).location, (3, 8))
        self.assertEquals(f.fragment_chunk(chunks[2]).location, (9, 6 + len(self.root_sequence)))

    def test_annotations(self):
        f = self.root.update('Bar')
        f.insert_bases(3, 'gataca')
        chunks = [chunk for chunk in f.chunks()]
        f.annotate(3, 8, 'A1', 'gene', 1)
        annotations = f.fragment_chunk(chunks[1]).annotations()
        self.assertEquals(len(annotations), 1)
        self.assertEquals(annotations[0].base_first, 3)
        self.assertEquals(annotations[0].base_last, 8)
        self.assertEquals(annotations[0].feature.name, 'A1')
        self.assertEquals(annotations[0].feature_base_first, 1)
        self.assertEquals(annotations[0].feature_base_last, 6)
