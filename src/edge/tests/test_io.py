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
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_fasta(f.name)
            filehandle = open(f.name, 'r')
            fasta = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        self.assertEquals(fasta, '>Bar\nagttcgaggctga\n')

    def test_outputs_gff_correctly(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', 1)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tgene\t2\t9\t.\t+\t.\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_can_output_gff_feature_on_reverse_strand(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', -1)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tgene\t2\t9\t.\t-\t.\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_after_insertion_correctly(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', 1)
        fragment.insert_bases(3, 'gataca')

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 19
Bar\tfeature\tgene\t2\t2\t.\t+\t.\tname=A1%%5B1:1%%5D
Bar\tfeature\tgene\t9\t15\t.\t+\t.\tname=A1%%5B2:8%%5D
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_with_hierarchical_annotations(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'gene', 1)
        fragment.annotate(2, 9, 'A1', 'mRNA', 1)
        fragment.annotate(2, 9, 'A1', 'CDS', 1)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tgene\t2\t9\t.\t+\t.\tname=A1
Bar\tfeature\tmRNA\t2\t9\t.\t+\t.\tname=A1
Bar\tfeature\tCDS\t2\t9\t.\t+\t0\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_features_with_null_phase(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'CDS', 1, qualifiers=dict(phase=None))

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tCDS\t2\t9\t.\t+\t0\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_features_with_integer_phase(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'CDS', 1, qualifiers=dict(phase=2))

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tCDS\t2\t9\t.\t+\t2\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_gff_features_with_phase_stored_as_array(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'CDS', 1, qualifiers=dict(phase=["2"]))

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tCDS\t2\t9\t.\t+\t2\tname=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)

    def test_outputs_qualifier_data_in_gff_features(self):
        fragment = self.fragment.indexed_fragment()
        fragment.annotate(2, 9, 'A1', 'CDS', 1,
                          qualifiers=dict(locus_tag="b0002",
                                          foo=["blah"],
                                          aliases=["foo", "bar"]))

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.close()
            IO(self.genome).to_gff(f.name)
            filehandle = open(f.name, 'r')
            gff = filehandle.read()
            filehandle.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region Bar 1 13
Bar\tfeature\tCDS\t2\t9\t.\t+\t0\taliases=foo,bar;foo=blah;locus_tag=b0002;name=A1
##FASTA
>Bar <unknown description>
%s
""" % (fragment.sequence,)
        self.assertEquals(expected, gff)
