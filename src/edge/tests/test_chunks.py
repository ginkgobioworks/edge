from unittest import TestCase
from edge import Connector, Fragment, FragmentNotFound
from edge.models import Edge


class FragmentTests(TestCase):

    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')
        self.root_sequence = 'agttcgaggctga'
        self.root = self.world.create_fragment_with_sequence('Foo', self.root_sequence)

    def test_root_fragment(self):
        self.assertEquals(self.root.sequence, self.root_sequence)
        self.assertEquals(self.root.name, 'Foo')

    def test_split_in_middle(self):
        u = self.root.update('Bar')
        self.assertEquals(len([c for c in u.chunks()]), 1)
        prev, cur = u._find_and_split_before(4)
        self.assertEquals([c.sequence for c in u.chunks()], ['agt', 'tcgaggctga'])
        self.assertEquals(u.sequence, self.root_sequence)

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
        prev, cur = u._find_and_split_before(len(self.root_sequence)+1)
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

    def test_get_sequence_by_region(self):
        for i in range(0, len(self.root_sequence)):
            for j in range(i, len(self.root_sequence)):
                self.assertEquals(self.root.get_sequence(i+1, j+1), self.root_sequence[i:j+1])

    def test_get_sequence_by_region_after_split(self):
        u = self.root.update('Bar')
        prev, cur = u._find_and_split_before(6)
        f = u.commit()
        for i in range(0, len(self.root_sequence)):
            for j in range(i, len(self.root_sequence)):
                self.assertEquals(f.get_sequence(i+1, j+1), self.root_sequence[i:j+1])

    def test_insert_sequence_in_middle(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gataca'+self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_start(self):
        u = self.root.update('Bar')
        u.insert_bases(1, 'gataca')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'gataca'+self.root_sequence)
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_end(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence+'gataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u.commit()
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence+'gatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence+'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_does_not_allow_insert_same_fragment_twice_which_creates_loops(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.insert_fragment(3, new_f)
        self.assertRaises(Exception, u.insert_fragment, 3, new_f)

    def test_insert_new_fragment_in_middle(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.insert_fragment(3, new_f)
        f = u.commit()
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gcccataca'+self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_start(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.insert_fragment(1, new_f)
        f = u.commit()
        self.assertEquals(f.sequence, 'gcccataca'+self.root_sequence)
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_end(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.insert_fragment(None, new_f)
        f = u.commit()
        self.assertEquals(f.sequence, self.root_sequence+'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_new_fragment_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.insert_fragment(None, new_f)
        f1 = u.commit()
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence+'gcccatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence+'gcccataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_existing_fragment_in_middle(self):
        existing_f = self.world.create_fragment_with_sequence('Test', 'gataca')
        u = self.root.update('Bar')
        u.insert_fragment(3, existing_f)
        f = u.commit()
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gataca'+self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_fragment_inherits_annotations_from_new_fragment(self):
        u = self.root.update('Bar')
        self.assertEquals([a for a in u.annotations() if a.name == 'Uma'], [])
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.annotate()
        new_f.annotate(2, 4, 'Uma', 'feature', 1)
        new_f = new_f.commit()
        u.insert_fragment(2, new_f)
        f = u.commit()
        self.assertEquals([a for a in f.annotations() if a.name == 'Uma'][0].first_bp, 3)
        self.assertEquals([a for a in f.annotations() if a.name == 'Uma'][0].last_bp, 5)

    def test_remove_sequence_in_middle(self):
        u = self.root.update('Bar')
        u.remove_bases(3, 4)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_start(self):
        u = self.root.update('Bar')
        u.remove_bases(1, 4)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[4:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_end(self):
        u = self.root.update('Bar')
        # removing before 4th to last bp, remove 4 bases
        u.remove_bases(len(self.root_sequence)-3, 4)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_past_end(self):
        u = self.root.update('Bar')
        # removing before 2nd to last bp, remove 4 bases
        u.remove_bases(len(self.root_sequence)-1, 4)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-2])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        # removing before 4th to last bp, remove 4 bases
        u.remove_bases(len(self.root_sequence)-3, 4)
        f1 = u.commit()
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence[:-4]+'aaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-4])
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_in_middle(self):
        u = self.root.update('Bar')
        u.replace_bases(3, 6, 'cccc')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'cccc'+self.root_sequence[8:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_start(self):
        u = self.root.update('Bar')
        u.replace_bases(1, 6, 'cccc')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'cccc'+self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_end(self):
        u = self.root.update('Bar')
        # removing before 6th to last bp, remove 6 bases
        u.replace_bases(len(self.root_sequence)-5, 6, 'cccc')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-6]+'cccc')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_past_end(self):
        u = self.root.update('Bar')
        # removing before 4th to last bp, remove 6 bases
        u.replace_bases(len(self.root_sequence)-3, 6, 'cccc')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4]+'cccc')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_sequence_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        # removing before 6th to last bp, remove 6 bases
        u.replace_bases(len(self.root_sequence)-5, 6, 'cccc')
        f1 = u.commit()
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence[:-6]+'ccccaaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-6]+'cccc')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_in_middle(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.replace_with_fragment(3, 6, new_f)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gcccataca'+self.root_sequence[8:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_start(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.replace_with_fragment(1, 6, new_f)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'gcccataca'+self.root_sequence[6:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_end(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.replace_with_fragment(len(self.root_sequence)-5, 6, new_f)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-6]+'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_past_end(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.replace_with_fragment(len(self.root_sequence)-3, 6, new_f)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[:-4]+'gcccataca')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_replace_fragment_at_end_then_insert_again(self):
        u = self.root.update('Bar')
        new_f = u._connector.create_fragment_with_sequence('Test', 'gataca')
        new_f = new_f.update('Test')
        new_f.insert_bases(2, 'ccc')
        new_f = new_f.commit()
        u.replace_with_fragment(len(self.root_sequence)-5, 6, new_f)
        f1 = u.commit()
        u = f1.update('Bar')
        u.insert_bases(None, 'aaaa')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence[:-6]+'gcccatacaaaaa')
        self.assertEquals(f1.sequence, self.root_sequence[:-6]+'gcccataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_multiple_operations(self):
        u = self.root.update('Bar')
        # removing before 4th to last bp, remove 6 bases
        u.replace_bases(len(self.root_sequence)-3, 6, 'cccc')
        u.insert_bases(None, 'ggg')
        u.remove_bases(1, 3)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        s = self.root_sequence
        self.assertEquals(f.sequence, s[3:-4]+'cccc'+'ggg')
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_uncommit(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        self.assertRaises(FragmentNotFound, Fragment, self.world, u.id)
        u.commit()
        f = Fragment(self.world, u.id)
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gataca'+self.root_sequence[2:])

    def test_insert_on_chunk_size_of_one(self):
        u = self.root.update('Bar')
        # this op creates a chunk with size of 1, then G, then the rest of the sequence
        u.insert_bases(2, 'G')
        # this op has to insert on the chunk with size of 1
        u.insert_bases(1, 'C')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, 'C'+self.root_sequence[0:1]+'G'+self.root_sequence[1:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_double_insert(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'G')
        u.insert_bases(3, 'C')
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'CG'+self.root_sequence[2:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_double_remove(self):
        u = self.root.update('Bar')
        u.remove_bases(3, 4)
        u.remove_bases(3, 4)
        f = u.commit()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.sequence, self.root_sequence[0:2]+self.root_sequence[10:])
        # does not affect root
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_insert_at_end_for_two_fragments(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u.commit()
        u = f1.update('Baz')
        u.insert_bases(None, 'atatat')
        f2 = u.commit()
        self.assertEquals(f2.sequence, self.root_sequence+'gataca'+'atatat')
        self.assertEquals(f1.sequence, self.root_sequence+'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_remove_then_insert(self):
        u = self.root.update('Bar')
        u.insert_bases(None, 'gataca')
        f1 = u.commit()
        u = f1.update('Baz')
        u.remove_bases(len(self.root_sequence)+1, 6)
        f2 = u.commit()
        u = f2.update('Far')
        u.insert_bases(None, 'atatat')
        f3 = u.commit()
        self.assertEquals(f3.sequence, self.root_sequence+'atatat')
        self.assertEquals(f2.sequence, self.root_sequence)
        self.assertEquals(f1.sequence, self.root_sequence+'gataca')
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_updater_indexes_chunk_locations_on_commit(self):
        from sqlalchemy.sql import select
        from edge.models import fragment_chunk_location_table

        # has one chunk before
        stmt = select([fragment_chunk_location_table.c.fragment_id])
        r = [x for x in self.world.conn.execute(stmt)]
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, self.root.id)

        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()

        # has three new chunks now
        stmt = select([fragment_chunk_location_table.c.fragment_id,
                       fragment_chunk_location_table.c.chunk_id,
                       fragment_chunk_location_table.c.base_first,
                       fragment_chunk_location_table.c.base_last])
        r = [x for x in self.world.conn.execute(stmt)]

        # each new chunk has an location index, additionally, the original chunk 1
        # has been divided into two chunks, each with an index.

        self.assertEquals(len(r), 5)
        self.assertItemsEqual(r, [
            (self.root.id, 1, 1, 2),
            (self.root.id, 2, 3, len(self.root_sequence)),
            (f.id, 1, 1, 2),
            (f.id, 3, 3, 8),
            (f.id, 2, 9, 6+len(self.root_sequence)),
        ])

    def test_add_edges_replaces_existing_edge_with_same_fragment_id(self):
        u = self.root.update('Bar')
        prev, cur = u._find_and_split_before(4)
        prev_chunk = u.get_chunk(prev)
        self.assertEquals(len(prev_chunk.out_edges), 1)
        self.assertEquals(prev_chunk.out_edges[0].to_chunk_id, cur)
        self.assertEquals(prev_chunk.out_edges[0].fragment_id, self.root.id)

        # add new edge, does not overwrite old edge since new edge has a new fragment id
        u._add_edges(prev, Edge(prev, 101, 100))
        prev_chunk = u.get_chunk(prev)
        self.assertEquals(len(prev_chunk.out_edges), 2)
        self.assertEquals(prev_chunk.out_edges[0].to_chunk_id, cur)
        self.assertEquals(prev_chunk.out_edges[0].fragment_id, self.root.id)
        self.assertEquals(prev_chunk.out_edges[1].to_chunk_id, 100)
        self.assertEquals(prev_chunk.out_edges[1].fragment_id, 101)

        # add another edge with fragment id 101, does overwrite existing edge with same fragment id
        u._add_edges(prev, Edge(prev, 101, 102))
        prev_chunk = u.get_chunk(prev)
        self.assertEquals(len(prev_chunk.out_edges), 2)
        self.assertEquals(prev_chunk.out_edges[0].to_chunk_id, cur)
        self.assertEquals(prev_chunk.out_edges[0].fragment_id, self.root.id)
        self.assertEquals(prev_chunk.out_edges[1].to_chunk_id, 102)
        self.assertEquals(prev_chunk.out_edges[1].fragment_id, 101)

    def test_find_chunk_by_location_index(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        f = f.annotate()

        # there should now be 3 chunks for f
        self.assertEquals([c.sequence for c in f.chunks()],
                          [self.root_sequence[0:2], 'gataca', self.root_sequence[2:]])
        chunk_ids = [c.id for c in f.chunks()]

        # get first bp
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_location_index(1)
        self.assertEquals(prev_chunk_id, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk_id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get last bp in first chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_location_index(2)
        self.assertEquals(prev_chunk_id, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk_id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get first bp in second chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_location_index(3)
        self.assertEquals(prev_chunk_id, chunk_ids[0])
        self.assertEquals(chunk.id, chunk_ids[1])
        self.assertEquals(next_chunk_id, chunk_ids[2])
        self.assertEquals(bases_visited, 2+6)

        # get last bp in last chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_location_index(6+len(self.root_sequence))
        self.assertEquals(prev_chunk_id, chunk_ids[1])
        self.assertEquals(chunk.id, chunk_ids[2])
        self.assertEquals(next_chunk_id, None)
        self.assertEquals(bases_visited, 6+len(self.root_sequence))

        # get bp past last chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_location_index(100)
        self.assertEquals(prev_chunk_id, chunk_ids[2])
        self.assertEquals(chunk, None)
        self.assertEquals(next_chunk_id, None)
        self.assertEquals(bases_visited, 6+len(self.root_sequence))

    def test_find_chunk_by_walking(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        f = f.annotate()

        # there should now be 3 chunks for f
        self.assertEquals([c.sequence for c in f.chunks()],
                          [self.root_sequence[0:2], 'gataca', self.root_sequence[2:]])
        chunk_ids = [c.id for c in f.chunks()]

        # get first bp
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_walking(1)
        self.assertEquals(prev_chunk_id, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk_id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get last bp in first chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_walking(2)
        self.assertEquals(prev_chunk_id, None)
        self.assertEquals(chunk.id, chunk_ids[0])
        self.assertEquals(next_chunk_id, chunk_ids[1])
        self.assertEquals(bases_visited, 2)

        # get first bp in second chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_walking(3)
        self.assertEquals(prev_chunk_id, chunk_ids[0])
        self.assertEquals(chunk.id, chunk_ids[1])
        self.assertEquals(next_chunk_id, chunk_ids[2])
        self.assertEquals(bases_visited, 2+6)

        # get last bp in last chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_walking(6+len(self.root_sequence))
        self.assertEquals(prev_chunk_id, chunk_ids[1])
        self.assertEquals(chunk.id, chunk_ids[2])
        self.assertEquals(next_chunk_id, None)
        self.assertEquals(bases_visited, 6+len(self.root_sequence))

        # get bp past last chunk
        prev_chunk_id, chunk, next_chunk_id, bases_visited = \
            f._find_chunk_prev_next_by_walking(100)
        self.assertEquals(prev_chunk_id, chunk_ids[2])
        self.assertEquals(chunk, None)
        self.assertEquals(next_chunk_id, None)
        self.assertEquals(bases_visited, 6+len(self.root_sequence))


class FragmentContextTest(TestCase):

    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')
        self.root_sequence = 'agttcgaggctga'
        self.root = self.world.create_fragment_with_sequence('Foo', self.root_sequence)

    def test_commits_after_with_block_and_can_access_child_fragment_via_last_updated(self):
        self.assertEquals(self.root.last_updated(), None)
        with self.root.update('Bar') as u:
            u.insert_bases(3, 'gataca')
        f = self.root.last_updated()
        self.assertEquals(f.name, 'Bar')
        self.assertEquals(f.parent_id, self.root.id)
        self.assertEquals(f.sequence, self.root_sequence[0:2]+'gataca'+self.root_sequence[2:])
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_does_not_commit_if_exception_in_block(self):
        class MyTestException(Exception):
            pass
        self.assertEquals(self.root.last_updated(), None)
        try:
            with self.root.update('Bar') as u:
                u.insert_bases(3, 'gataca')
                raise MyTestException('error')
        except MyTestException:
            pass
        self.assertEquals(self.root.last_updated(), None)
        self.assertEquals(self.root.sequence, self.root_sequence)

    def test_commits_but_does_not_create_updated_fragment_when_annotating(self):
        self.assertEquals(self.root.last_updated(), None)
        self.assertEquals(len(self.root.annotations()), 0)
        with self.root.annotate() as u:
            u.annotate(3, 4, 'A', 'gene', 1)
        self.assertEquals(self.root.last_updated(), None)
        self.assertEquals(self.root.sequence, self.root_sequence)
        self.assertEquals(len(self.root.annotations()), 1)


class FragmentChunkTest(TestCase):

    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')
        self.root_sequence = 'agttcgaggctga'
        self.root = self.world.create_fragment_with_sequence('Foo', self.root_sequence)

    def test_next_chunk_id(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(chunks[0].next_chunk_id, chunks[1].id)
        self.assertEquals(chunks[1].next_chunk_id, chunks[2].id)
        self.assertEquals(chunks[2].next_chunk_id, None)

    def test_next_chunk(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(chunks[0].next_chunk.id, chunks[1].id)
        self.assertEquals(chunks[1].next_chunk.id, chunks[2].id)
        self.assertEquals(chunks[2].next_chunk, None)

    def test_prev_chunk(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(chunks[0].prev_chunk, None)
        self.assertEquals(chunks[1].prev_chunk.id, chunks[0].id)
        self.assertEquals(chunks[2].prev_chunk.id, chunks[1].id)

    def test_location(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        chunks = [chunk for chunk in f.chunks()]
        self.assertEquals(chunks[0].location, (1, 2))
        self.assertEquals(chunks[1].location, (3, 8))
        print chunks[2].sequence
        self.assertEquals(chunks[2].location, (9, 6+len(self.root_sequence)))

    def test_annotations(self):
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        f = u.commit()
        chunks = [chunk for chunk in f.chunks()]
        a = f.annotate()
        a.annotate(3, 8, 'A1', 'gene', 1)
        f = a.commit()
        annotations = chunks[1].annotations()
        self.assertEquals(len(annotations), 1)
        self.assertEquals(annotations[0].first_bp, None)
        self.assertEquals(annotations[0].last_bp, None)
        self.assertEquals(annotations[0].name, 'A1')
        self.assertEquals(annotations[0].annotation_first_bp, 1)
        self.assertEquals(annotations[0].annotation_last_bp, 6)
