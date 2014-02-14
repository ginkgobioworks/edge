from django.test import TestCase
from edge.models import *
from edge.io import IO
import tempfile
import os


class IOTest(TestCase):

    def setUp(self):
        self.root_sequence = 'agttcgaggctga'
        self.root = Fragment.create_with_sequence('Foo', self.root_sequence)

    def test_outputs_fasta(self):
        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(self.root).to_fasta(f.name)
            h = open(f.name, 'r')
            fasta = h.read()
            h.close()
            os.unlink(f.name)

        self.assertEquals(fasta, '>fragment %s|Foo\nagttcgaggctga\n' % (self.root.id,))

    def test_outputs_gff(self):
        f = self.root.indexed_fragment()
        f.annotate(2, 9, 'A1', 'gene', 1)
        u = f.update('Bar')
        u.insert_bases(3, 'gataca')
        frag = u

        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(frag).to_gff(f.name)
            h = open(f.name, 'r')
            gff = h.read()
            print gff
            h.close()
            os.unlink(f.name)

        # be aware of the tabs in the string below
        expected = """##gff-version 3
##sequence-region fragment %(fid)s: Bar 1 19
fragment %(fid)s: Bar\tfeature\tgene\t2\t2\t.\t+\t.\tname=A1
fragment %(fid)s: Bar\tfeature\tgene\t9\t15\t.\t+\t.\tname=A1
""" % {'fid': u.id}
        self.assertEquals(expected, gff)
