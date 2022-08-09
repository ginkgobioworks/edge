from django.test import TestCase

from edge.models import Genome


class GenomeTest(TestCase):
    def test_create_genome(self):
        self.assertEquals(len(Genome.objects.all()), 0)
        genome = Genome.create("Foo")
        self.assertEquals(len(Genome.objects.all()), 1)
        self.assertEquals(Genome.objects.all()[0].name, "Foo")
        self.assertEquals(Genome.objects.all()[0].id, genome.id)

    def test_add_fragments_to_genome_in_place(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        self.assertEquals(len(genome.fragments.all()), 1)
        self.assertEquals(genome.fragments.all()[0].name, "chrI")
        self.assertEquals(genome.fragments.all()[0].indexed_fragment().sequence, s)

    def test_add_fragments_to_genome_not_in_place_creates_and_updates_child(self):
        parent = Genome.create("Foo")
        self.assertEquals(len(Genome.objects.all()), 1)
        u = parent.update()
        s = "atggcatattcgcagct"
        u.add_fragment("chrI", s)

        child = u
        # created a child genome
        self.assertCountEqual(
            [g.id for g in Genome.objects.all()], [parent.id, child.id]
        )
        # child genome points to parent genome
        self.assertEquals(child.parent.id, parent.id)
        # child genome has changes, parent does not
        self.assertEquals(len(child.fragments.all()), 1)
        self.assertEquals(len(parent.fragments.all()), 0)

    def test_find_annotation_by_name(self):
        genome = Genome.create("Foo")
        s = "atggcatattcgcagct"
        f = genome.add_fragment("chrI", s)
        f.annotate(3, 8, "Foo gene", "gene", 1)

        annotations = genome.indexed_genome().find_annotation_by_name("Foo gene")
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

        # case insensitive
        annotations = genome.indexed_genome().find_annotation_by_name("foo gene")
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

        annotations = genome.indexed_genome().find_annotation_by_name("Foo bar")
        self.assertEquals(len(annotations), 0)

    def test_find_annotation_by_qualifier(self):
        genome = Genome.create("Foo")
        s = "atggcatattcgcagct"
        f = genome.add_fragment("chrI", s)
        f.annotate(
            3, 8, "Foo gene", "gene", 1, qualifiers=dict(foo="bar,baz", far="foo,fu")
        )

        # finds bar
        annotations = genome.indexed_genome().find_annotation_by_qualifier("bar")
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

        # case insensitive
        annotations = genome.indexed_genome().find_annotation_by_qualifier("bAr")
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

        # does not find bar,b
        annotations = genome.indexed_genome().find_annotation_by_qualifier("bar,b")
        self.assertEquals(len(annotations), 0)

        # can limit to a field
        annotations = genome.indexed_genome().find_annotation_by_qualifier(
            "bar", fields=["foo"]
        )
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

        # bar in qualifier, but not in 'far' field
        annotations = genome.indexed_genome().find_annotation_by_qualifier(
            "bar", fields=["far"]
        )
        self.assertEquals(len(annotations), 0)

        # can search in multiple fields
        annotations = genome.indexed_genome().find_annotation_by_qualifier(
            "bar", fields=["far", "foo"]
        )
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

    def test_find_annotation_by_feature(self):
        genome = Genome.create("Foo")
        s = "atggcatattcgcagct"
        f = genome.add_fragment("chrI", s)
        feature = f.annotate(3, 8, "Foo gene", "gene", 1)

        annotations = genome.indexed_genome().find_annotation_by_feature(feature)
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].base_first, 3)
        self.assertEquals(annotations[f.id][0].base_last, 8)
        self.assertEquals(annotations[f.id][0].feature.name, "Foo gene")

    def test_changes_return_empty_array_if_no_parent(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        self.assertEquals(genome.indexed_genome().changes(), [])

    def test_can_insert_and_get_changes(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)

        # insert
        g_u = genome.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual(
            [c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)]
        )

    def test_update_fragment_by_id(self):
        genome = Genome.create("Foo")
        s = "atggcatattcgcagct"
        f1 = genome.add_fragment("chrI", s)

        g2 = genome.update()
        with g2.update_fragment_by_fragment_id(f1.id) as f2:
            f2.insert_bases(3, "gataca")

        self.assertEquals(
            g2.fragments.all()[0].indexed_fragment().sequence, s[0:2] + "gataca" + s[2:]
        )
        self.assertEquals(genome.fragments.all()[0].indexed_fragment().sequence, s)

    def test_can_update_fragment_by_name_and_assign_new_name(self):
        genome = Genome.create("Foo")
        genome.add_fragment("chrI", "atggcatattcgcagct")
        self.assertCountEqual([f.name for f in genome.fragments.all()], ["chrI"])

        # insert
        g = genome.update()
        with g.update_fragment_by_name("chrI", "foobar") as f:
            f.insert_bases(3, "gataca")

        self.assertCountEqual([fr.name for fr in g.fragments.all()], ["foobar"])

    def test_can_update_fragment_by_id_and_assign_new_name(self):
        genome = Genome.create("Foo")
        f0 = genome.add_fragment("chrI", "atggcatattcgcagct")
        self.assertCountEqual([f.name for f in genome.fragments.all()], ["chrI"])

        # insert
        g = genome.update()
        with g.update_fragment_by_fragment_id(f0.id, "foobar") as f:
            f.insert_bases(3, "gataca")

        self.assertCountEqual([fr.name for fr in g.fragments.all()], ["foobar"])

    def test_can_insert_then_insert_and_get_second_insert_only_as_changes(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        g1 = genome

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
        g2 = g_u

        # insert again
        g_u = g2.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(9, "gataca")
        g3 = g_u

        changes = g2.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual(
            [c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)]
        )

        changes = g3.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual(
            [c.location for c in changes], [(3, 8), (9, 14), (15, len(s) + 6 + 6)]
        )

    def test_only_get_changes_from_changed_fragment(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        genome.add_fragment("chrII", s)
        g1 = genome

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
        g2 = g_u

        changes = g2.indexed_genome().changes()
        self.assertEquals(len(changes), 3)
        self.assertCountEqual(
            [c.location for c in changes], [(1, 2), (3, 8), (9, len(s) + 6)]
        )

    def test_can_remove_and_get_changes_back(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        g = genome

        # remove
        g_u = g.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.remove_bases(3, 4)
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 2)
        self.assertCountEqual([c.location for c in changes], [(1, 2), (3, len(s) - 4)])

    def test_can_insert_and_remove_and_get_all_changes_back(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        g = genome

        # insert and remove
        g_u = g.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
            f.remove_bases(10, 4)
        g = g_u

        changes = g.indexed_genome().changes()
        self.assertEquals(len(changes), 4)
        self.assertCountEqual(
            [c.location for c in changes],
            [(1, 2), (3, 8), (9, 9), (10, len(s) + 6 - 4)],
        )

    def test_can_add_notes_on_update(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)

        # insert
        u = genome.update(notes=u"Bar bar")
        with u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
        u = Genome.objects.get(pk=u.pk)
        self.assertEquals(u.notes, u"Bar bar")

    def test_can_set_name_on_update(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)

        # insert
        u = genome.update(name=u"Bar bar")
        with u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
        u = Genome.objects.get(pk=u.pk)
        self.assertEquals(u.name, u"Bar bar")

    def test_get_changed_locations_by_fragment(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        s0 = "atggcatattcgcagct"
        s1 = "gtacggctagtcgatt"
        s2 = "acgatcgggattgagtcgattc"

        # add initial sequence
        f1 = genome.add_fragment("chrI", s0 + s1 + s2)
        f2 = genome.add_fragment("chrII", s0 + s1 + s2)

        # annotate it to break it up into chunks
        with genome.annotate_fragment_by_name("chrI") as f:
            f.annotate(1, len(s0), "F1", "feature", 1)
            f.annotate(len(s0) + 1, len(s0) + len(s1), "F2", "feature", 1)
            f.annotate(
                len(s0) + len(s1) + 1, len(s0) + len(s1) + len(s2), "F3", "feature", 1
            )

        with genome.annotate_fragment_by_name("chrII") as f:
            f.annotate(1, len(s0), "F1", "feature", 1)
            f.annotate(len(s0) + 1, len(s0) + len(s1), "F2", "feature", 1)
            f.annotate(
                len(s0) + len(s1) + 1, len(s0) + len(s1) + len(s2), "F3", "feature", 1
            )

        # insert
        u = genome.update()
        with u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
            f1 = f
        with u.update_fragment_by_name("chrII") as f:
            f.insert_bases(3, "gataca")
            f.insert_bases(6 + len(s0) + len(s1) + 2, "gataca")
            f2 = f

        g2 = Genome.objects.get(pk=u.pk)
        changes = g2.indexed_genome().changed_locations_by_fragment()
        for f in changes:
            if f.id == f1.id:
                self.assertEquals(changes[f], [[1, len(s0) + 6]])
            elif f.id == f2.id:
                self.assertEquals(
                    changes[f],
                    [
                        [1, len(s0) + 6],
                        [
                            len(s0) + len(s1) + 6 + 1,
                            len(s0) + 6 + len(s1) + len(s2) + 6,
                        ],
                    ],
                )
            else:
                raise Exception("Unexpected fragment")

    def test_get_coordinate_diff_from_parent_genome(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        g1 = genome.indexed_genome()
        f1_id = g1.fragments.first().id

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(3, "gataca")
            f.remove_bases(10, 4)
        g2 = g_u.indexed_genome()
        f2_id = g2.fragments.first().id

        # insert again
        g_u = g2.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(9, "gataca")
        g3 = g_u.indexed_genome()
        f3_id = g3.fragments.first().id

        # insert again inside of previous integration
        g_u = g3.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(5, "gataca")
        g4 = g_u.indexed_genome()
        f4_id = g4.fragments.first().id

        diff1 = g2.get_coordinate_diff_from_parent_genome(g1)
        truediff1 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 3, 'parent_ends_before': 3,
             'child_fragment_name': 'chrI', 'child_fragment_id': f2_id,
             'child_starts_at': 3, 'child_ends_before': 9},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 4, 'parent_ends_before': 8,
             'child_fragment_name': 'chrI', 'child_fragment_id': f2_id,
             'child_starts_at': 10, 'child_ends_before': 10}
        ]
        self.assertEqual(diff1, truediff1)

        diff2 = g3.get_coordinate_diff_from_parent_genome(g2)
        truediff2 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f2_id,
             'parent_starts_at': 9, 'parent_ends_before': 9,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 9, 'child_ends_before': 15}
        ]
        self.assertEqual(diff2, truediff2)

        diff3 = g4.get_coordinate_diff_from_parent_genome(g3)
        truediff3 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f3_id,
             'parent_starts_at': 5, 'parent_ends_before': 5,
             'child_fragment_name': 'chrI', 'child_fragment_id': f4_id,
             'child_starts_at': 5, 'child_ends_before': 11}
        ]
        self.assertEqual(diff3, truediff3)

        diff4 = g4.get_coordinate_diff_from_parent_genome(g2)
        truediff4 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f2_id,
             'parent_starts_at': 5, 'parent_ends_before': 5,
             'child_fragment_name': 'chrI', 'child_fragment_id': f4_id,
             'child_starts_at': 5, 'child_ends_before': 11},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f2_id,
             'parent_starts_at': 9, 'parent_ends_before': 9,
             'child_fragment_name': 'chrI', 'child_fragment_id': f4_id,
             'child_starts_at': 15, 'child_ends_before': 21}
        ]
        self.assertEqual(diff4, truediff4)

        diff5 = g3.get_coordinate_diff_from_parent_genome(g1)
        truediff5 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 3, 'parent_ends_before': 3,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 3, 'child_ends_before': 15},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 4, 'parent_ends_before': 8,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 16, 'child_ends_before': 16}
        ]
        self.assertEqual(diff5, truediff5)

        diff6 = g4.get_coordinate_diff_from_parent_genome(g1)
        truediff6 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 3, 'parent_ends_before': 3,
             'child_fragment_name': 'chrI', 'child_fragment_id': f4_id,
             'child_starts_at': 3, 'child_ends_before': 21},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 4, 'parent_ends_before': 8,
             'child_fragment_name': 'chrI', 'child_fragment_id': f4_id,
             'child_starts_at': 22, 'child_ends_before': 22},
        ]
        self.assertEqual(diff6, truediff6)

    def test_get_coordinate_diff_from_parent_genome_at_beginning_and_end(self):
        genome = Genome.create("Foo")
        self.assertEquals(len(genome.fragments.all()), 0)

        # add initial sequence
        s = "atggcatattcgcagct"
        genome.add_fragment("chrI", s)
        g1 = genome.indexed_genome()
        f1_id = g1.fragments.first().id

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.insert_bases(17, "atatat")
            f.insert_bases(10, "atcg")
            f.insert_bases(1, "gcgcgc")
        g2 = g_u.indexed_genome()
        f2_id = g2.fragments.first().id

        # insert
        g_u = g2.update()
        with g_u.update_fragment_by_name("chrI") as f:
            f.remove_bases(27, 6)
            f.remove_bases(1, 6)
        g3 = g_u.indexed_genome()
        f3_id = g3.fragments.first().id

        diff1 = g2.get_coordinate_diff_from_parent_genome(g1)
        truediff1 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 1, 'parent_ends_before': 1,
             'child_fragment_name': 'chrI', 'child_fragment_id': f2_id,
             'child_starts_at': 1, 'child_ends_before': 7},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 10, 'parent_ends_before': 10,
             'child_fragment_name': 'chrI', 'child_fragment_id': f2_id,
             'child_starts_at': 16, 'child_ends_before': 20},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 17, 'parent_ends_before': 17,
             'child_fragment_name': 'chrI', 'child_fragment_id': f2_id,
             'child_starts_at': 27, 'child_ends_before': 33},
        ]
        self.assertEqual(diff1, truediff1)

        diff2 = g3.get_coordinate_diff_from_parent_genome(g2)
        truediff2 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f2_id,
             'parent_starts_at': 1, 'parent_ends_before': 7,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 1, 'child_ends_before': 1},
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f2_id,
             'parent_starts_at': 27, 'parent_ends_before': 33,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 21, 'child_ends_before': 21},
        ]
        self.assertEqual(diff2, truediff2)

        diff3 = g3.get_coordinate_diff_from_parent_genome(g1)
        truediff3 = [
            {'parent_fragment_name': 'chrI', 'parent_fragment_id': f1_id,
             'parent_starts_at': 10, 'parent_ends_before': 10,
             'child_fragment_name': 'chrI', 'child_fragment_id': f3_id,
             'child_starts_at': 10, 'child_ends_before': 14},
        ]
        self.assertEqual(diff3, truediff3)
