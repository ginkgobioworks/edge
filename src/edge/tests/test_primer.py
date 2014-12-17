from django.test import TestCase
from edge.primer import primer3_design
from edge.models import Genome, Fragment, Genome_Fragment


class Primer3Test(TestCase):

    def build_genome(self, circular, *templates):
        g = Genome(name='Foo')
        g.save()
        for seq in templates:
            f = Fragment.create_with_sequence('Bar', seq, circular=circular)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
        return g

    def test_primer3_finds_primers(self):
        upstream = "cagtacgatcgttggtatgctgactactagcgtagctagcacgtcgtgtccaggcttgagcgacgt"
        product = "cagctggtaatcgtactcgtactagcatcgtacgtgtctgatcatctgacgtatcatctga"
        downstream = "agtgacgtcgtgtgtagcgtactgtatcgtgtgtcgcgcgtagtcatctgatcgtacgtactgaat"
        template = ''.join([upstream, product, downstream])

        print 'up %s, len %s' % (len(upstream), len(product))
        g = self.build_genome(False, 'a'*40+template+'a'*40)
        f = g.fragments.all()[0]

        r = primer3_design(f, 40+len(upstream)+1, len(product), 50, 50, {})
        print r
        self.assertEquals(True, False)
