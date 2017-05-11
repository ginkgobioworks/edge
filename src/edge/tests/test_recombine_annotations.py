import os

from Bio.Seq import Seq
from django.test import TestCase

import edge.orfs
from edge.recombine import find_swap_region_with_annotations, recombine
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
        return Genome.objects.get(pk=g.id)

    def setUp(self):
        self.upstream = "gagattgtccgcgtttt"
        self.front_bs = "catagcgcacaggacgcggag"
        self.middle = "cggcaccttaattgcgaattgcgagctgacgtctgcatgtagccg"
        self.back_bs = "taatgaccccgaagcagg"
        self.downstream = "gttaaggcgcgaacat"
        self.template = ''.join([self.upstream, self.front_bs, self.middle,
                                 self.back_bs, self.downstream])
        self.arm_len = min(len(self.front_bs), len(self.back_bs))
        self.genome = self.build_genome(False, self.template)
        self.fragment = self.genome.fragments.all()[0].indexed_fragment()

        self.old_min_protein_len = edge.orfs.min_protein_len
        edge.orfs.min_protein_len = 10

    def tearDown(self):
        edge.orfs.min_protein_len = self.old_min_protein_len

    def test_does_not_preserve_annotation_if_replaced_is_different(self):
        replaced = "aaaaaaaaaaaaaaaaaaa"
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs),
                               len(self.upstream) + len(self.front_bs) + len(self.middle) - 1,
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 0)

    def test_preserves_annotation_with_single_mutation(self):
        replaced = [c for c in self.middle]
        self.assertEquals(replaced[10], 'a')
        replaced[10] = 'c'
        replaced = ''.join(replaced)
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle))
        self.assertEquals(a['feature_name'], 'Foo A11C')

    def test_detects_insertion(self):
        replaced = self.middle[0:10] + 'c' + self.middle[10:]
        self.assertEquals(replaced[10], 'c')
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle) + 1)
        self.assertEquals(a['feature_name'], 'Foo +11C')

    def test_detects_deletion(self):
        replaced = self.middle[0:13] + self.middle[14:]
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle) - 1)
        self.assertEquals(a['feature_name'], 'Foo -14G')

    def test_preserves_unchanged_annotation(self):
        replaced = self.middle[0:10] + self.middle[11:]
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 12,
                               len(self.upstream) + len(self.front_bs) + len(self.middle) - 1,
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 12 - 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle) - 1 - 1)
        self.assertEquals(a['feature_name'], 'Foo')

    def test_preserves_feature_direction_and_type(self):
        replaced = self.middle[0:13] + self.middle[14:]
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'foobar', -1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle) - 1)
        self.assertEquals(a['feature_name'], 'Foo -14G')
        self.assertEquals(a['feature_type'], 'foobar')
        self.assertEquals(a['feature_strand'], -1)

    def test_preserves_multiple_annotations(self):
        replaced = self.middle[0:13] + self.middle[14:]
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + 2, len(self.upstream) + 10, 'Bar', 'static', -1)
        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'changed', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 2)

        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], 2)
        self.assertEquals(a['base_last'], 10)
        self.assertEquals(a['feature_name'], 'Bar')
        self.assertEquals(a['feature_type'], 'static')
        self.assertEquals(a['feature_strand'], -1)

        a = r[0].cassette_annotations[1]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle) - 1)
        self.assertEquals(a['feature_name'], 'Foo -14G')
        self.assertEquals(a['feature_type'], 'changed')
        self.assertEquals(a['feature_strand'], 1)

    def test_computes_annotation_for_reverse_cassette(self):
        replaced = [c for c in self.middle]
        self.assertEquals(replaced[10], 'a')
        replaced[10] = 'c'
        replaced = ''.join(replaced)
        cassette = ''.join([self.front_bs, replaced, self.back_bs])
        cassette = str(Seq(cassette).reverse_complement())

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)
        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(self.middle))
        self.assertEquals(a['feature_name'], 'Foo A11C')

    def test_returns_new_orf(self):
        replaced = 'atgatcatcatcatcatcatcatcatcatcatcatcatcatcatcatctag'
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)

        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(replaced))
        self.assertEquals(a['feature_name'], 'ORF frame 1')
        self.assertEquals(a['feature_type'], 'ORF')
        self.assertEquals(a['feature_strand'], 1)

    def test_finds_orf_in_reverse_direction(self):
        replaced = 'atgatcatcatcatcatcatcatcatcatcatcatcatcatcatcatctag'
        replaced = str(Seq(replaced).reverse_complement())
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        r = find_swap_region_with_annotations(self.genome, cassette, self.arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].cassette_annotations), 1)

        a = r[0].cassette_annotations[0]
        self.assertEquals(a['base_first'], len(self.front_bs) + 1)
        self.assertEquals(a['base_last'], len(self.front_bs) + len(replaced))
        self.assertEquals(a['feature_name'], 'ORF frame 1')
        self.assertEquals(a['feature_type'], 'ORF')
        self.assertEquals(a['feature_strand'], -1)

    def test_adds_multiple_annotations_to_modified_genome(self):
        replaced = self.middle[0:13] + self.middle[14:]
        cassette = ''.join([self.front_bs, replaced, self.back_bs])

        self.fragment.annotate(len(self.upstream) + 2, len(self.upstream) + 10, 'Bar', 'static', -1)
        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'changed', 1)

        c = recombine(self.genome, cassette, self.arm_len)
        f = c.fragments.all()[0].indexed_fragment()

        annotations = f.annotations()
        self.assertEquals(len(annotations), 3)

        a = annotations[0]
        self.assertEquals(a.feature.type, 'operation')

        a = annotations[1]
        self.assertEquals(a.feature.name, 'Bar')
        self.assertEquals(a.feature.type, 'static')
        self.assertEquals(a.feature.strand, -1)
        self.assertEquals(a.base_first, len(self.upstream) + 2)
        self.assertEquals(a.base_last, len(self.upstream) + 10)

        a = annotations[2]
        self.assertEquals(a.feature.name, 'Foo -14G')
        self.assertEquals(a.feature.type, 'changed')
        self.assertEquals(a.feature.strand, 1)
        self.assertEquals(a.base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(a.base_last, len(
            self.upstream) + len(self.front_bs) + len(self.middle) - 1)

    def test_adds_annotations_correctly_with_reverse_cassette(self):
        replaced = [c for c in self.middle]
        self.assertEquals(replaced[10], 'a')
        replaced[10] = 'c'
        replaced = ''.join(replaced)
        cassette = ''.join([self.front_bs, replaced, self.back_bs])
        cassette = str(Seq(cassette).reverse_complement())

        self.fragment.annotate(len(self.upstream) + len(self.front_bs) + 1,
                               len(self.upstream) + len(self.front_bs) + len(self.middle),
                               'Foo', 'gene', 1)

        c = recombine(self.genome, cassette, self.arm_len)
        f = c.fragments.all()[0].indexed_fragment()

        annotations = f.annotations()
        self.assertEquals(len(annotations), 2)

        a = annotations[0]
        self.assertEquals(a.feature.type, 'operation')

        a = annotations[1]
        self.assertEquals(a.feature.name, 'Foo A11C')
        self.assertEquals(a.feature.type, 'gene')
        self.assertEquals(a.feature.strand, 1)
        self.assertEquals(a.base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(a.base_last, len(self.upstream) + len(self.front_bs) + len(self.middle))
