from unittest import TestCase
from edge import Connector, Genome
import os
import tempfile


class GenomeTest(TestCase):

    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')

    def test_create_genome(self):
        self.assertEquals(len(self.world.genomes()), 0)
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(self.world.genomes()), 1)
        self.assertEquals(self.world.genomes()[0].name, 'Foo')
        self.assertEquals(self.world.genomes()[0].id, genome.id)

    def test_non_genomic_fragments(self):
        genome = self.world.create_genome('Foo')
        u = genome._edit()
        s = 'atggcatattcgcagct'
        f1 = u.add_fragment('chrI', s)
        u.commit()
        f2 = self.world.create_fragment_with_sequence('Bar', 'aacctaaaattataa')
        self.assertEquals(len(self.world.non_genomic_fragments()), 1)
        self.assertEquals(self.world.non_genomic_fragments()[0].name, 'Bar')
        self.assertEquals(self.world.non_genomic_fragments()[0].id, f2.id)

    def test_add_fragments_to_genome_in_place(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)
        u = genome._edit()
        s = 'atggcatattcgcagct'
        f = u.add_fragment('chrI', s)
        u.commit()
        self.assertEquals(len(genome.fragments()), 1)
        self.assertEquals(genome.fragments()[0].name, 'chrI')
        self.assertEquals(genome.fragments()[0].sequence, s)

    def test_genome_updates_not_visible_after_fragment_commit_and_before_genome_commit(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)
        u = genome._edit()
        s = 'atggcatattcgcagct'
        f = u.add_fragment('chrI', s)
        f.commit()
        # have not yet committed genome updater
        self.assertEquals(len(genome.fragments()), 0)
        # commit genome updater
        u.commit()
        self.assertEquals(len(genome.fragments()), 1)

    def test_add_fragments_to_genome_not_in_place_creates_and_updates_child(self):
        parent = self.world.create_genome('Foo')
        self.assertEquals(len(self.world.genomes()), 1)
        u = parent.update()
        s = 'atggcatattcgcagct'
        f = u.add_fragment('chrI', s)
        child = u.commit()
        # refetch from db
        child = Genome(self.world, child.id)

        # created a child genome
        self.assertItemsEqual([g.id for g in self.world.genomes()], [parent.id, child.id])
        # child genome points to parent genome
        self.assertEquals(child.parent.id, parent.id)
        # child genome has changes, parent does not
        self.assertEquals(len(child.fragments()), 1)
        self.assertEquals(len(parent.fragments()), 0)

    def test_import_gff_creates_fragments_and_annotate_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t30\t80\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
chrII\tTest\tgene\t40\t60\t.\t-\t.\tID=f4;gene=g4
chrII\tTest\tgene\t20\t80\t.\t+\t.\tID=i5;Name=f5
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
>chrII
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        genome = self.world.create_genome('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()

            u = genome._edit()
            u.import_gff(f.name)
            u.commit()

            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertItemsEqual([f.name for f in genome.fragments()], ['chrI', 'chrII'])

        # verify chrI fragment
        chrI = [f for f in genome.fragments() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].first_bp, 20)
        self.assertEquals(chrI.annotations()[1].last_bp, 28)
        self.assertEquals(chrI.annotations()[1].name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].strand, 1)
        self.assertEquals(chrI.annotations()[0].first_bp, 30)
        self.assertEquals(chrI.annotations()[0].last_bp, 80)
        self.assertEquals(chrI.annotations()[0].name, 'f2')
        self.assertEquals(chrI.annotations()[0].strand, -1)

        # verify chrII fragment
        chrII = [f for f in genome.fragments() if f.name == 'chrII'][0]
        self.assertEquals(len(chrII.sequence), 160)
        # consecutive annotations merged even though they span multiple chunks
        self.assertEquals(len(chrII.annotations()), 2)
        self.assertEquals(chrII.annotations()[0].first_bp, 40)
        self.assertEquals(chrII.annotations()[0].last_bp, 60)
        self.assertEquals(chrII.annotations()[0].name, 'g4')  # has gene, use gene name
        self.assertEquals(chrII.annotations()[0].strand, -1)
        self.assertEquals(chrII.annotations()[1].first_bp, 20)
        self.assertEquals(chrII.annotations()[1].last_bp, 80)
        self.assertEquals(chrII.annotations()[1].name, 'f5')
        self.assertEquals(chrII.annotations()[1].strand, 1)

    def test_find_annotation(self):
        genome = self.world.create_genome('Foo')
        u = genome._edit()
        s = 'atggcatattcgcagct'
        f = u.add_fragment('chrI', s)
        f.annotate(3, 8, 'Foo gene', 'gene', 1)
        u.commit()
        annotations = genome.find_annotation('Foo gene')
        self.assertEquals(len(annotations), 1)
        self.assertEquals(f.id in annotations, True)
        self.assertEquals(len(annotations[f.id]), 1)
        self.assertEquals(annotations[f.id][0].first_bp, 3)
        self.assertEquals(annotations[f.id][0].last_bp, 8)
        self.assertEquals(annotations[f.id][0].name, 'Foo gene')

    def test_changes_return_empty_array_if_no_parent(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)
        u = genome._edit()
        s = 'atggcatattcgcagct'
        f = u.add_fragment('chrI', s)
        g = u.commit()
        self.assertEquals(g.changes(), [])

    def test_can_insert_and_get_changes(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        g_u = genome._edit()
        s = 'atggcatattcgcagct'
        g_u.add_fragment('chrI', s)
        g = g_u.commit()

        # insert
        g_u = g.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g = g_u.commit()

        changes = g.changes()
        self.assertEquals(len(changes), 3)
        self.assertItemsEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s)+6)])

    def test_can_insert_then_insert_and_get_second_insert_only_as_changes(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        g_u = genome._edit()
        s = 'atggcatattcgcagct'
        g_u.add_fragment('chrI', s)
        g1 = g_u.commit()

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g2 = g_u.commit()

        # insert again
        g_u = g2.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(9, 'gataca')
        g3 = g_u.commit()

        changes = g2.changes()
        self.assertEquals(len(changes), 3)
        self.assertItemsEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s)+6)])

        changes = g3.changes()
        self.assertEquals(len(changes), 3)
        self.assertItemsEqual([c.location for c in changes], [(3, 8), (9, 14), (15, len(s)+6+6)])

    def test_only_get_changes_from_changed_fragment(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        g_u = genome._edit()
        s = 'atggcatattcgcagct'
        g_u.add_fragment('chrI', s)
        g_u.add_fragment('chrII', s)
        g1 = g_u.commit()

        # insert
        g_u = g1.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
        g2 = g_u.commit()

        changes = g2.changes()
        self.assertEquals(len(changes), 3)
        self.assertItemsEqual([c.location for c in changes], [(1, 2), (3, 8), (9, len(s)+6)])

    def test_can_remove_and_get_changes_back(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        g_u = genome._edit()
        s = 'atggcatattcgcagct'
        g_u.add_fragment('chrI', s)
        g = g_u.commit()

        # remove
        g_u = g.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.remove_bases(3, 4)
        g = g_u.commit()

        changes = g.changes()
        self.assertEquals(len(changes), 2)
        self.assertItemsEqual([c.location for c in changes], [(1, 2), (3, len(s)-4)])

    def test_can_insert_and_remove_and_get_all_changes_back(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        g_u = genome._edit()
        s = 'atggcatattcgcagct'
        g_u.add_fragment('chrI', s)
        g = g_u.commit()

        # insert and remove
        g_u = g.update()
        with g_u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')
            f.remove_bases(10, 4)
        g = g_u.commit()

        changes = g.changes()
        self.assertEquals(len(changes), 4)
        self.assertItemsEqual([c.location for c in changes],
                              [(1, 2), (3, 8), (9, 9), (10, len(s)+6-4)])

    def test_can_add_notes_on_update(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        with genome._edit() as u:
            s = 'atggcatattcgcagct'
            u.add_fragment('chrI', s)

        # insert
        with genome.update(notes=u'Bar bar') as u:
            with u.update_fragment_by_name('chrI') as f:
                f.insert_bases(3, 'gataca')
        g = genome.last_updated()
        self.assertEquals(g.notes, u'Bar bar')

    def test_can_set_name_on_update(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        # add initial sequence
        with genome._edit() as u:
            s = 'atggcatattcgcagct'
            u.add_fragment('chrI', s)

        # insert
        with genome.update(name=u'Bar bar') as u:
            with u.update_fragment_by_name('chrI') as f:
                f.insert_bases(3, 'gataca')
        g = genome.last_updated()
        self.assertEquals(g.name, u'Bar bar')

    def test_get_changed_locations_by_fragment(self):
        genome = self.world.create_genome('Foo')
        self.assertEquals(len(genome.fragments()), 0)

        s0 = 'atggcatattcgcagct'
        s1 = 'gtacggctagtcgatt'
        s2 = 'acgatcgggattgagtcgattc'

        # add initial sequence
        with genome._edit() as u:
            u.add_fragment('chrI', s0+s1+s2)
            u.add_fragment('chrII', s0+s1+s2)
        # annotate it to break it up into chunks
        with genome.annotate() as u:
            with u.annotate_fragment_by_name('chrI') as f:
                f.annotate(1, len(s0), 'F1', 'feature', 1)
                f.annotate(len(s0)+1, len(s0)+len(s1), 'F2', 'feature', 1)
                f.annotate(len(s0)+len(s1)+1, len(s0)+len(s1)+len(s2), 'F3', 'feature', 1)
            with u.annotate_fragment_by_name('chrII') as f:
                f.annotate(1, len(s0), 'F1', 'feature', 1)
                f.annotate(len(s0)+1, len(s0)+len(s1), 'F2', 'feature', 1)
                f.annotate(len(s0)+len(s1)+1, len(s0)+len(s1)+len(s2), 'F3', 'feature', 1)

        # insert
        with genome.update() as u:
            with u.update_fragment_by_name('chrI') as f:
                f.insert_bases(3, 'gataca')
            with u.update_fragment_by_name('chrII') as f:
                f.insert_bases(3, 'gataca')
                f.insert_bases(6+len(s0)+len(s1)+2, 'gataca')

        g2 = genome.last_updated()
        changes = g2.changed_locations_by_fragment()
        for f in changes:
            if f.id == 3:
                self.assertEquals(changes[f], [[1, len(s0)+6]])
            elif f.id == 4:
                self.assertEquals(changes[f],
                                  [[1, len(s0)+6],
                                   [len(s0)+len(s1)+6+1, len(s0)+6+len(s1)+len(s2)+6]])
            else:
                raise Exception('Unexpected fragment')


class GenomeContextTest(TestCase):
    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')
        self.genome = self.world.create_genome('Foo')

    def test_edit_in_place_does_not_return_updated_genome(self):
        self.assertEquals(self.genome.last_updated(), None)
        with self.genome._edit() as u:
            s = 'atggcatattcgcagct'
            f = u.add_fragment('chrI', s)
        self.assertEquals(self.genome.last_updated(), None)
        self.assertEquals(len(self.genome.fragments()), 1)
        self.assertEquals(self.genome.fragments()[0].name, 'chrI')
        self.assertEquals(self.genome.fragments()[0].sequence, s)

    def test_edit_in_place_does_not_commit_if_exception_raised(self):
        class MyTestException(Exception):
            pass
        try:
            with self.genome._edit() as u:
                s = 'atggcatattcgcagct'
                f = u.add_fragment('chrI', s)
                raise MyTestException('error')
        except MyTestException:
            pass
        self.assertEquals(self.genome.last_updated(), None)
        self.assertEquals(len(self.genome.fragments()), 0)

    def test_cannot_update_fragment_if_editing_or_annotating_genome(self):
        with self.genome._edit() as u:
            s = 'atggcatattcgcagct'
            f = u.add_fragment('chrI', s)
        with self.genome._edit() as u:
            try:
                with u.update_fragment_by_name('chrI') as f:
                    raise Exception('Should not be here')
            except:
                pass
        with self.genome.annotate() as u:
            try:
                with u.update_fragment_by_name('chrI') as f:
                    raise Exception('Should not be here')
            except:
                pass
        with self.genome.update() as u:
            with u.update_fragment_by_name('chrI') as f:
                # this works!
                pass

    def test_update_commits_and_can_get_updated_genome(self):
        self.assertEquals(self.genome.last_updated(), None)
        with self.genome.update() as u:
            s = 'atggcatattcgcagct'
            f = u.add_fragment('chrI', s)
        self.assertNotEqual(self.genome.last_updated(), None)
        g = self.genome.last_updated()
        self.assertEquals(len(g.fragments()), 1)
        self.assertEquals(g.fragments()[0].name, 'chrI')
        self.assertEquals(g.fragments()[0].sequence, s)

    def test_update_does_not_commit_if_exception_raised(self):
        self.assertEquals(self.genome.last_updated(), None)

        class MyTestException(Exception):
            pass

        try:
            with self.genome.update() as u:
                s = 'atggcatattcgcagct'
                f = u.add_fragment('chrI', s)
                raise MyTestException('error')
        except MyTestException:
            pass
        self.assertEquals(self.genome.last_updated(), None)
        self.assertEquals(len(self.genome.fragments()), 0)

    def test_update_fragment_by_name(self):
        self.assertEquals(self.genome.last_updated(), None)
        with self.genome.update() as u:
            s = 'atggcatattcgcagct'
            f = u.add_fragment('chrI', s)
        g1 = self.genome.last_updated()
        with g1.update() as u:
            with u.update_fragment_by_name('chrI') as f:
                f.insert_bases(3, 'gataca')
        g2 = g1.last_updated()

        self.assertEquals(g2.fragments()[0].sequence, s[0:2]+'gataca'+s[2:])
        self.assertEquals(g1.fragments()[0].sequence, s)

    def test_update_fragment_by_id(self):
        self.assertEquals(self.genome.last_updated(), None)
        with self.genome.update() as u:
            s = 'atggcatattcgcagct'
            f1 = u.add_fragment('chrI', s)
        g1 = self.genome.last_updated()
        with g1.update() as u:
            with u.update_fragment_by_fragment_id(f1.id) as f2:
                f2.insert_bases(3, 'gataca')
        g2 = g1.last_updated()

        self.assertEquals(g2.fragments()[0].sequence, s[0:2]+'gataca'+s[2:])
        self.assertEquals(g1.fragments()[0].sequence, s)
