import os
from django.test import TestCase
from edge.primer import design_primers
from edge.models import Genome, Fragment, Genome_Fragment
from edge.pcr import pcr_from_genome
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn


class DesignPrimerTest(TestCase):

    def build_genome(self, circular, *templates):
        g = Genome(name='Foo')
        g.save()
        for seq in templates:
            f = Fragment.create_with_sequence('Bar', seq, circular=circular)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
            try:
                os.unlink(fragment_fasta_fn(f))
            except:
                pass
        build_all_genome_dbs(refresh=True)
        return g

    def test_primer3_finds_primers_when_given_range_to_find_primers(self):
        upstream = "cagtacgatcgttggtatgctgactactagcgtagctagcacgtcgtgtccaggcttgagcgacgt"
        product = "cagctggtaatcgtactcgtactagcatcgtacgtgtctgatcatctgacgtatcatctga"
        downstream = "agtgacgtcgtgtgtagcgtactgtatcgtgtgtcgcgcgtagtcatctgatcgtacgtactgaat"
        template = ''.join([upstream, product, downstream])

        g = self.build_genome(False, 'a'*40+template+'a'*40)
        f = g.fragments.all()[0]

        res = design_primers(f, 40+len(upstream)+1, len(product), 50, 50, {})
        self.assertEquals(len(res), 5)

        for r in res:
            p = pcr_from_genome(g, r['PRIMER_LEFT_SEQUENCE'], r['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(product) >= 0, True)
