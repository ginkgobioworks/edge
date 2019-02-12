import os
import tempfile

from django.test import TestCase

from edge.io import IO
from edge.models import (
    Fragment,
    Genome,
    Genome_Fragment,
)


class IOTest(TestCase):

    def setUp(self):
        self.sequence = 'agttcgaggctga'
        self.genome = Genome.create('Foo')
        self.fragment = Fragment.create_with_sequence('Bar', self.sequence)
        Genome_Fragment(genome=self.genome, fragment=self.fragment, inherited=False).save()

    def test_outputs_fasta(self):
        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(self.genome).to_fasta(f.name)
            h = open(f.name, 'r')
            fasta = h.read()
            h.close()
            os.unlink(f.name)

        self.assertEquals(fasta, '>Bar\nagttcgaggctga\n')

    def test_outputs_gff(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', 1)
        fragment.insert_bases(3, 'gataca')

        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            h = open(f.name, 'r')
            gff = h.read()
            h.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 19
Bar\tfeature\tgene\t2\t2\t.\t+\t.\tname=A1
Bar\tfeature\tgene\t9\t15\t.\t+\t.\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_handles_negative_strand(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', -1)
        fragment.insert_bases(3, 'gataca')

        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            h = open(f.name, 'r')
            gff = h.read()
            h.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 19
Bar\tfeature\tgene\t2\t2\t.\t-\t.\tname=A1
Bar\tfeature\tgene\t9\t15\t.\t-\t.\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)
