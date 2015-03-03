import os
import json
from django.test import TestCase
from Bio.Seq import Seq
from edge.recombine import find_swap_region
from edge.models import Genome, Fragment, Genome_Fragment
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn


class GenomeRecombinationAnnotationsTest(TestCase):
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

    def setUp(self):
        self.upstream = "gagattgtccgcgtttt"
        self.front_bs = "catagcgcacaggacgcggag"
        self.middle = "cggcaccttaattgcgagattgcgagctgacgtctgcatgtgagccg"
        self.back_bs = "taatgaccccgaagcagg"
        self.downstream = "gttaaggcgcgaacat"
        self.template = ''.join([self.upstream, self.front_bs, self.middle,
                                 self.back_bs, self.downstream])
        self.arm_len = min(len(self.front_bs), len(self.back_bs))
        self.genome = self.build_genome(False, self.template)
        self.fragment = self.genome.fragments.all()[0].indexed_fragment()

    def test_does_not_preserve_annotation_if_replaced_is_different(self):
        replaced = "aaaaaaaaaaaaaaaaaaa"
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream)+len(self.front_bs),
                               len(self.upstream)+len(self.front_bs)+len(self.middle)-1,
                               'Foo', 'gene', 1)

        r = find_swap_region(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 0)

    def test_preserves_annotation_with_single_mutation(self):
        replaced = [c for c in self.middle]
        self.assertEquals(replaced[10], 'a')
        replaced[10] = 'c'
        replaced = ''.join(replaced)

        cassette = ''.join([self.front_bs, replaced, self.back_bs])
        self.fragment.annotate(len(self.upstream)+len(self.front_bs),
                               len(self.upstream)+len(self.front_bs)+len(self.middle)-1,
                               'Foo', 'gene', 1)

        r = find_swap_region(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs))
        self.assertEquals(a['base_last'], len(self.front_bs)+len(self.middle)-1)
        self.assertEquals(a['feature_name'], 'Foo A11C')

    def test_detects_insertion(self):
        self.assertEquals(True, False)

    def test_detects_deletion(self):
        self.assertEquals(True, False)

    def test_computes_annotation_for_reverse_cassette(self):
        self.assertEquals(True, False)
    
    def test_preserves_unchanged_annotation(self):
        self.assertEquals(True, False)

    def test_preserves_multiple_annotations(self):
        self.assertEquals(True, False)
    
    def test_returns_new_orf(self):
        self.assertEquals(True, False)

    def test_adds_multiple_annotations_to_modified_genome(self):
        self.assertEquals(True, False)
