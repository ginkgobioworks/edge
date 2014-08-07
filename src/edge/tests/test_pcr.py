import os
import json
from Bio.Seq import Seq
from django.test import TestCase
from edge.op import pcr_from_genome


class GenomePcrTest(TestCase):

    def test_pcr_produces_expected_product(self):
        from edge.models import Genome, Fragment, Genome_Fragment
        from edge.management.commands.build_edge_blastdb import build_fragment_db, fragment_fasta_fn

        upstream = "gagattgtccgcgtttt"
        primer_1_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        primer_2_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        template = ''.join([upstream, primer_1_bs, middle, primer_2_bs, downstream])

        primer_1 = 'aaaaaaaaaa'+primer_1_bs
        primer_2 = 'tttttttttt'+'cctgcttcggggtcatta'

        g = Genome(name='Foo')
        g.save()
        f = Fragment.create_with_sequence('Bar', template)
        Genome_Fragment(genome=g, fragment=f, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f))
        except:
            pass
        build_fragment_db()

        r = pcr_from_genome(g, primer_1, primer_2)
        self.assertEquals(r[0], ''.join([primer_1, middle, str(Seq(primer_2).reverse_complement())]))
