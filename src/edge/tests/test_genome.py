from django.test import TestCase

from edge.models import Genome


class GenomeTest(TestCase):

    def test_create_genome(self):
        self.assertEquals(len(Genome.objects.all()), 0)
        genome = Genome.create('Foo')
        self.assertEquals(len(Genome.objects.all()), 1)
        self.assertEquals(Genome.objects.all()[0].name, 'Foo')
        self.assertEquals(Genome.objects.all()[0].id, genome.id)

    def test_add_fragments_to_genome_in_place(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        self.assertEquals(len(genome.fragments.all()), 1)
        self.assertEquals(genome.fragments.all()[0].name, 'chrI')
        self.assertEquals(genome.fragments.all()[0].indexed_fragment().sequence, s)

    def test_add_fragments_to_genome_not_in_place_creates_and_updates_child(self):
        parent = Genome.create('Foo')
        self.assertEquals(len(Genome.objects.all()), 1)
        u = parent.update()
        s = 'atggcatattcgcagct'
        u.add_fragment('chrI', s)

        child = u
        # created a child genome
        self.assertCountEqual([g.id for g in Genome.objects.all()], [parent.id, child.id])
        # child genome points to parent genome
        self.assertEquals(child.parent.id, parent.id)
        # child genome has changes, parent does not
        self.assertEquals(len(child.fragments.all()), 1)
        self.assertEquals(len(parent.fragments.all()), 0)

    def test_find_annotation_by_name(self):
        genome = Genome.create('Foo')
        s = 'atggcatattcgcagct'
        f = genome.add_fragment('chrI', s)
        f.annotate(3, 8, 'Foo gene', 'gene', 1)

        annotations = genome.indexed_genome().find_annotation_by_name('Foo gene')
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

        # case insensitive
        annotations = genome.indexed_genome().find_annotation_by_name('foo gene')
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

        annotations = genome.indexed_genome().find_annotation_by_name('Foo bar')
        self.assertEquals(len(annotations), 0)

    def test_find_annotation_by_qualifier(self):
        genome = Genome.create('Foo')
        s = 'atggcatattcgcagct'
        f = genome.add_fragment('chrI', s)
        f.annotate(3, 8, 'Foo gene', 'gene', 1, qualifiers=dict(foo='bar,baz', far='foo,fu'))

        # finds bar
        annotations = genome.indexed_genome().find_annotation_by_qualifier('bar')
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

        # case insensitive
        annotations = genome.indexed_genome().find_annotation_by_qualifier('bAr')
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

        # does not find bar,b
        annotations = genome.indexed_genome().find_annotation_by_qualifier('bar,b')
        self.assertEquals(len(annotations), 0)

        # can limit to a field
        annotations = genome.indexed_genome().find_annotation_by_qualifier('bar', fields=['foo'])
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

        # bar in qualifier, but not in 'far' field
        annotations = genome.indexed_genome().find_annotation_by_qualifier('bar', fields=['far'])
        self.assertEquals(len(annotations), 0)

        # can search in multiple fields
        annotations = \
            genome.indexed_genome().find_annotation_by_qualifier('bar', fields=['far', 'foo'])
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

    def test_find_annotation_by_feature(self):
        genome = Genome.create('Foo')
        s = 'atggcatattcgcagct'
        f = genome.add_fragment('chrI', s)
        feature = f.annotate(3, 8, 'Foo gene', 'gene', 1)

        annotations = genome.indexed_genome().find_annotation_by_feature(feature)
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, 'Foo gene')

    def test_changes_return_empty_array_if_no_parent(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        self.assertEquals(genome.indexed_genome().changes(), [])

    def test_can_insert_and_get_changes(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)

        # insert
        g_u = genome.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)])

    def test_update_fragment_by_id(self):
        genome = Genome.create('Foo')
        s = 'atggcatattcgcagct'
        f1 = genome.add_fragment('chrI', s)

        g2 = genome.update()
        with g2.update_fragment_by_fragment_id(f1.id) as f2:
            f2.insert_bases(3, 'gataca')

        self.assertEquals(g2.fragments.all()[0].indexed_fragment(
        ).sequence, s[0:2] + 'gataca' + s[2:])
        self.assertEquals(genome.fragments.all()[0].indexed_fragment().sequence, s)

    def test_can_update_fragment_by_name_and_assign_new_name(self):
        genome = Genome.create('Foo')
        genome.add_fragment('chrI', 'atggcatattcgcagct')
        self.assertCountEqual([f.name for f in genome.fragments.all()], ['chrI'])

        # insert
        g = genome.update()
        with g.update_fragment_by_name('chrI', 'foobar') as f:
            f.insert_bases(3, 'gataca')

        self.assertCountEqual([fr.name for fr in g.fragments.all()], ['foobar'])

    def test_can_update_fragment_by_id_and_assign_new_name(self):
        genome = Genome.create('Foo')
        f0 = genome.add_fragment('chrI', 'atggcatattcgcagct')
        self.assertCountEqual([f.name for f in genome.fragments.all()], ['chrI'])

        # insert
        g = genome.update()
        with g.update_fragment_by_fragment_id(f0.id, 'foobar') as f:
            f.insert_bases(3, 'gataca')

        self.assertCountEqual([fr.name for fr in g.fragments.all()], ['foobar'])

    def test_can_insert_then_insert_and_get_second_insert_only_as_changes(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        g1 = genome

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g2 = g_u

        # insert again
        g_u = g2.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(9, 'gataca')
        g3 = g_u

        changes = g2.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)])

        changes = g3.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual([c.location for c in changes], [
                              (3, 8), (9, 14), (15, len(s) + 6 + 6)])

    def test_only_get_changes_from_changed_fragment(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        genome.add_fragment('chrII', s)
        g1 = genome

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g2 = g_u

        changes = g2.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)])

    def test_can_remove_and_get_changes_back(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        g = genome

        # remove
        g_u = g.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.remove_bases(3, 4)
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 2)
        self.assertCountEqual([c.location for c in changes], [(1, 2), (3, len(s) - 4)])

    def test_can_insert_and_remove_and_get_all_changes_back(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)
        g = genome

        # insert and remove
        g_u = g.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
            f.remove_bases(10, 4)
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 4)
        self.assertCountEqual([c.location for c in changes],
                              [(1, 2), (3, 8), (9, 9), (10, len(s) + 6 - 4)])

    def test_can_add_notes_on_update(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)

        # insert
        u = genome.update(notes=u'Bar bar')
        with u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        u = Genome.objects.get(pk=u.pk)
        self.assertEquals(u.notes, u'Bar bar')

    def test_can_set_name_on_update(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = 'atggcatattcgcagct'
        genome.add_fragment('chrI', s)

        # insert
        u = genome.update(name=u'Bar bar')
        with u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        u = Genome.objects.get(pk=u.pk)
        self.assertEquals(u.name, u'Bar bar')

    def test_get_changed_locations_by_fragment(self):
        genome = Genome.create('Foo')
        self.assertEquals(len(genome.fragments.all()), 0)

        s0 = 'atggcatattcgcagct'
        s1 = 'gtacggctagtcgatt'
        s2 = 'acgatcgggattgagtcgattc'

        # add initial sequence
        f1 = genome.add_fragment('chrI', s0 + s1 + s2)
        f2 = genome.add_fragment('chrII', s0 + s1 + s2)

        # annotate it to break it up into chunks
        with genome.annotate_fragment_by_name('chrI') as f:
            f.annotate(1, len(s0), 'F1', 'feature', 1)
            f.annotate(len(s0) + 1, len(s0) + len(s1), 'F2', 'feature', 1)
            f.annotate(len(s0) + len(s1) + 1, len(s0) + len(s1) + len(s2), 'F3', 'feature', 1)

        with genome.annotate_fragment_by_name('chrII') as f:
            f.annotate(1, len(s0), 'F1', 'feature', 1)
            f.annotate(len(s0) + 1, len(s0) + len(s1), 'F2', 'feature', 1)
            f.annotate(len(s0) + len(s1) + 1, len(s0) + len(s1) + len(s2), 'F3', 'feature', 1)

        # insert
        u = genome.update()
        with u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
            f1 = f
        with u.update_fragment_by_name('chrII') as f:
            f.insert_bases(3, 'gataca')
            f.insert_bases(6 + len(s0) + len(s1) + 2, 'gataca')
            f2 = f

        g2 = Genome.objects.get(pk=u.pk)
        changes = g2.indexed_genome().changed_locations_by_fragment()
        for f in changes:
            if f.id == f1.id:
                self.assertEquals(changes[f], [[1, len(s0) + 6]])
            elif f.id == f2.id:
                self.assertEquals(
                    changes[f],
                    [[1, len(s0) + 6],
                     [len(s0) + len(s1) + 6 + 1, len(s0) + 6 + len(s1) + len(s2) + 6]])
            else:
                raise Exception('Unexpected fragment')
