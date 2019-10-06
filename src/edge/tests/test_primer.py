import os
from Bio.Seq import Seq
from django.test import TestCase
from edge.primer import design_primers, design_primers_from_template
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
            except OSError:
                pass
        build_all_genome_dbs(refresh=True)
        return Genome.objects.get(pk=g.id)

    def test_primer3_finds_primers_when_given_range_to_find_primers(self):
        upstream = 'cagtacgatcgttggtatgctgactactagcgtagctagcacgtcgtgtccaggcttgagcgacgt'
        product = 'cagctggtaatcgtactcgtactagcatcgtacgtgtctgatcatctgacgtatcatctga'
        downstream = 'agtgacgtcgtgtgtagcgtactgtatcgtgtgtcgcgcgtagtcatctgatcgtacgtactgaat'
        template = ''.join([upstream, product, downstream])

        g = self.build_genome(False, 'a' * 40 + template + 'a' * 40)
        f = g.fragments.all()[0]

        res = design_primers(f, 40 + len(upstream) + 1, len(product), 50, 50, {})
        self.assertEquals(len(res), 5)

        for r in res:
            p = pcr_from_genome(g, r['PRIMER_LEFT_SEQUENCE'], r['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(product) >= 0, True)

    def test_computes_primer_distance_to_junction(self):
        upstream = 'cagtacgatcgttggtatgctgactactagcgtagctagcacgtcgtgtccaggcttgagcgacgt'
        product = 'cagctggtaatcgtactcgtactagcatcgtacgtgtctgatcatctgacgtatcatctga'
        downstream = 'agtgacgtcgtgtgtagcgtactgtatcgtgtgtcgcgcgtagtcatctgatcgtacgtactgaat'
        template = ''.join([upstream, product, downstream])

        WIN = 10
        junctions = [len(upstream) + 1, len(upstream) + len(product) - 1]
        res = design_primers_from_template(template, len(
            upstream) + 1 - WIN, len(product) + WIN * 2, junctions, {})
        self.assertEquals(len(res), 5)
        for r in res:
            left = r['PRIMER_LEFT_SEQUENCE']
            self.assertEquals(r['PRIMER_LEFT_SEQUENCE_DISTANCE_TO_JUNCTION'] >= WIN, True)
            self.assertEquals(template.lower().index(left.lower())
                              + len(left)
                              + r['PRIMER_LEFT_SEQUENCE_DISTANCE_TO_JUNCTION'],
                              junctions[0])
            right = str(Seq(r['PRIMER_RIGHT_SEQUENCE']).reverse_complement())
            self.assertEquals(r['PRIMER_RIGHT_SEQUENCE_DISTANCE_TO_JUNCTION'] >= WIN, True)
            self.assertEquals(template.lower().index(right.lower())
                              - r['PRIMER_RIGHT_SEQUENCE_DISTANCE_TO_JUNCTION'],
                              junctions[1])
