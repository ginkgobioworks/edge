from django.test import TestCase
from django.db import IntegrityError

from edge.management.commands.remove_fragment import remove_fragment
from edge.management.commands.remove_genome import remove_genome
from edge.models import (
    Genome,
    Fragment,
)


class RemoveTest(TestCase):

    def setUp(self):
        import os
        import tempfile
        from edge import import_gff

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
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()
            import_gff('TestGenome', f.name)
            self.genome = Genome.objects.get(name='TestGenome')
            os.unlink(f.name)
        self.assertCountEqual([fr.name for fr in self.genome.fragments.all()], ['chrI', 'chrII'])

    def test_remove_fragment_works(self):
        self.assertCountEqual([f.name for f in self.genome.fragments.all()], ['chrI', 'chrII'])
        self.assertCountEqual([f.name for f in Fragment.objects.all()], ['chrI', 'chrII'])

        chrI = [f for f in self.genome.fragments.all() if f.name == 'chrI'][0]
        remove_fragment(chrI.id)

        self.assertCountEqual([f.name for f in self.genome.fragments.all()], ['chrII'])
        self.assertCountEqual([f.name for f in Fragment.objects.all()], ['chrII'])

    def test_cannot_remove_fragment_if_it_has_derived_fragment(self):

        chrI = [f.indexed_fragment() for f in self.genome.fragments.all() if f.name == 'chrI'][0]

        u = chrI.update('Bar')
        u.insert_bases(3, 'gataca')

        self.assertCountEqual([f.name for f in self.genome.fragments.all()], ['chrI', 'chrII'])
        self.assertCountEqual([f.name for f in Fragment.objects.all()], ['chrI', 'chrII', 'Bar'])

        self.assertRaises(IntegrityError, remove_fragment, chrI.id)

        self.assertCountEqual([f.name for f in self.genome.fragments.all()], ['chrI', 'chrII'])
        self.assertCountEqual([f.name for f in Fragment.objects.all()], ['chrI', 'chrII', 'Bar'])

    def test_remove_genome_works(self):
        self.assertCountEqual([g.id for g in Genome.objects.all()], [self.genome.id])
        remove_genome(self.genome.id)
        self.assertEquals(list(Genome.objects.all()), [])

    def test_cannot_remove_genome_if_it_has_derived_genome(self):

        u = self.genome.update()
        with u.update_fragment_by_name('chrI') as f:
            f.insert_bases(3, 'gataca')

        self.assertCountEqual([g.id for g in Genome.objects.all()], [self.genome.id, u.id])
        self.assertRaises(IntegrityError, remove_genome, self.genome.id)
        self.assertCountEqual([g.id for g in Genome.objects.all()], [self.genome.id, u.id])

        # can remove derived genome then parent genome
        remove_genome(u.id)
        remove_genome(self.genome.id)
        self.assertEquals(list(Genome.objects.all()), [])
