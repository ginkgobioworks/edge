from django.test import TestCase
from edge.models import *
import os
import tempfile


class ImporterTest(TestCase):

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

        genome = Genome.create('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()

            u = genome.edit()
            u.import_gff(f.name)

            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertItemsEqual([f.name for f in genome.fragments.all()], ['chrI', 'chrII'])

        # verify chrI fragment
        chrI = [f for f in genome.fragments.all() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(chrI.annotations()[1].feature.name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 30)
        self.assertEquals(chrI.annotations()[0].base_last, 80)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

        # verify chrII fragment
        chrII = [f for f in genome.fragments.all() if f.name == 'chrII'][0]
        self.assertEquals(len(chrII.sequence), 160)
        # consecutive annotations merged even though they span multiple chunks
        self.assertEquals(len(chrII.annotations()), 2)
        self.assertEquals(chrII.annotations()[0].base_first, 40)
        self.assertEquals(chrII.annotations()[0].base_last, 60)
        self.assertEquals(chrII.annotations()[0].feature.name, 'g4')  # has gene, use gene name
        self.assertEquals(chrII.annotations()[0].feature.strand, -1)
        self.assertEquals(chrII.annotations()[1].base_first, 20)
        self.assertEquals(chrII.annotations()[1].base_last, 80)
        self.assertEquals(chrII.annotations()[1].feature.name, 'f5')
        self.assertEquals(chrII.annotations()[1].feature.strand, 1)

    def test_import_feature_starting_at_first_base(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t1\t80\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        genome = Genome.create('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()
            u = genome.edit()
            u.import_gff(f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [f for f in genome.fragments.all() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(chrI.annotations()[1].feature.name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 1)
        self.assertEquals(chrI.annotations()[0].base_last, 80)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

    def test_import_partially_overlapping_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t19\t21\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        genome = Genome.create('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()
            u = genome.edit()
            u.import_gff(f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [f for f in genome.fragments.all() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(chrI.annotations()[1].feature.name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 19)
        self.assertEquals(chrI.annotations()[0].base_last, 21)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

    def test_import_overlapping_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t20\t28\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        genome = Genome.create('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()
            u = genome.edit()
            u.import_gff(f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [f for f in genome.fragments.all() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(chrI.annotations()[1].feature.name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 20)
        self.assertEquals(chrI.annotations()[0].base_last, 28)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

    def test_import_feature_ending_at_last_base(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t20\t28\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t160\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        genome = Genome.create('Foo')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write(data)
            f.close()
            u = genome.edit()
            u.import_gff(f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [f for f in genome.fragments.all() if f.name == 'chrI'][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 160)
        self.assertEquals(chrI.annotations()[1].feature.name, 'i3')  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 20)
        self.assertEquals(chrI.annotations()[0].base_last, 28)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)
