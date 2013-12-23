from unittest import TestCase
from edge import Connector, IO
import tempfile
import os


class IOTest(TestCase):

    def setUp(self):
        self.world = Connector.create_db('/tmp/world.db')
        self.root_sequence = 'agttcgaggctga'
        self.root = self.world.create_fragment_with_sequence('Foo', self.root_sequence)

    def test_outputs_fasta(self):
        with tempfile.NamedTemporaryFile(mode='rw+', delete=False) as f:
            f.close()
            IO(self.root).to_fasta(f.name)
            h = open(f.name, 'r')
            fasta = h.read()
            print h.read()
            h.close()
            os.unlink(f.name)

        self.assertEquals(fasta, '>fragment 1|Foo\nagttcgaggctga\n')

    def test_outputs_gff(self):
        a = self.root.annotate()
        a.annotate(2, 9, 'A1', 'gene', 1)
        a.commit()
        u = self.root.update('Bar')
        u.insert_bases(3, 'gataca')
        frag = u.commit()

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
##sequence-region fragment 2: Bar 1 19
fragment 2: Bar\tfeature\tgene\t2\t2\t.\t+\t.\tname=A1
fragment 2: Bar\tfeature\tgene\t9\t15\t.\t+\t.\tname=A1
"""
        self.assertEquals(expected, gff)
