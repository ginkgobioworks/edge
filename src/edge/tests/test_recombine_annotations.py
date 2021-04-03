import os
import json

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
            except OSError:
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

    def test_preserves_annotations_on_homology_arm_fwd_that_have_not_changed(self):
        cassette = ''.join([self.front_bs, self.back_bs])

        # add annotaiton on upstream arm
        self.fragment.annotate(len(self.upstream) + 2,
                               len(self.upstream) + len(self.front_bs),
                               "Up arm", "feature", 1)

        c = recombine(self.genome, cassette, self.arm_len)
        f = c.fragments.all()[0].indexed_fragment()

        annotations = f.annotations()
        self.assertEquals(len(annotations), 2)

        a = annotations[0]
        self.assertEquals(a.feature.type, 'operation')

        a = annotations[1]
        self.assertEquals(a.feature.name, 'Up arm')
        self.assertEquals(a.feature.type, 'feature')
        self.assertEquals(a.feature.strand, 1)
        self.assertEquals(a.base_first, len(self.upstream) + 2)
        self.assertEquals(a.base_last, len(self.upstream) + len(self.front_bs))

    def test_preserves_annotations_on_homology_arm_rev_that_have_not_changed(self):
        cassette = ''.join([self.front_bs, self.back_bs])

        # add annotaiton on reverse strand of downstream arm
        self.fragment.annotate(
            len(self.upstream) + len(self.front_bs) + len(self.middle) + 1,
            len(self.upstream) + len(self.front_bs) + len(self.middle) + len(self.back_bs),
            "Down arm", "feature", -1)

        c = recombine(self.genome, cassette, self.arm_len)
        f = c.fragments.all()[0].indexed_fragment()

        annotations = f.annotations()
        self.assertEquals(len(annotations), 2)

        a = annotations[0]
        self.assertEquals(a.feature.type, 'operation')

        a = annotations[1]
        self.assertEquals(a.feature.name, 'Down arm')
        self.assertEquals(a.feature.type, 'feature')
        self.assertEquals(a.feature.strand, -1)
        self.assertEquals(a.base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(a.base_last, len(self.upstream) + len(self.front_bs) + len(self.back_bs))

    def test_adds_annotations_on_fwd_strand(self):
        donor = "a" * 100 + "g" * 100 + "c" * 50 + "t" * 50 + "g" * 100 + "c" * 100
        cassette = ''.join([self.front_bs, donor, self.back_bs])
        flen = len(self.front_bs)

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        c = recombine(self.genome, cassette, self.arm_len, annotations=annotations)
        f = c.fragments.all()[0].indexed_fragment()
        fragment_sequence = f.sequence

        ans = f.annotations()
        self.assertEquals(len(ans), 7)

        ans = sorted(ans, key=lambda a: (a.base_first, -a.base_last))

        self.assertEquals(ans[0].feature.type, 'operation')

        self.assertEquals(ans[1].feature.type, 'promoter')
        self.assertEquals(ans[1].feature.name, 'pFavorite')
        self.assertEquals(ans[1].feature.strand, 1)
        self.assertEquals(ans[1].base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(ans[1].base_last, len(self.upstream) + len(self.front_bs) + 100)
        self.assertEquals(fragment_sequence[ans[1].base_first - 1:ans[1].base_last], "a" * 100)

        self.assertEquals(ans[2].feature.type, 'gene')
        self.assertEquals(ans[2].feature.name, 'Favorite')
        self.assertEquals(ans[2].feature.strand, 1)
        self.assertEquals(ans[2].base_first, len(self.upstream) + len(self.front_bs) + 100 + 1)
        self.assertEquals(ans[2].base_last, len(self.upstream) + len(self.front_bs) + 100 + 100)
        self.assertEquals(fragment_sequence[ans[2].base_first - 1:ans[2].base_last], "g" * 100)

        self.assertEquals(ans[3].feature.type, 'terminator')
        self.assertEquals(ans[3].feature.name, 'tFavorite')
        self.assertEquals(ans[3].feature.strand, 1)
        self.assertEquals(ans[3].base_first, len(self.upstream) + len(self.front_bs) + 200 + 1)
        self.assertEquals(ans[3].base_last, len(self.upstream) + len(self.front_bs) + 200 + 50)
        self.assertEquals(fragment_sequence[ans[3].base_first - 1:ans[3].base_last], "c" * 50)

        self.assertEquals(ans[4].feature.type, 'terminator')
        self.assertEquals(ans[4].feature.name, 'tBest')
        self.assertEquals(ans[4].feature.strand, -1)
        self.assertEquals(ans[4].base_first, len(self.upstream) + len(self.front_bs) + 250 + 1)
        self.assertEquals(ans[4].base_last, len(self.upstream) + len(self.front_bs) + 250 + 50)
        self.assertEquals(fragment_sequence[ans[4].base_first - 1:ans[4].base_last], "t" * 50)

        self.assertEquals(ans[5].feature.type, 'gene')
        self.assertEquals(ans[5].feature.name, 'Best')
        self.assertEquals(ans[5].feature.strand, -1)
        self.assertEquals(ans[5].base_first, len(self.upstream) + len(self.front_bs) + 300 + 1)
        self.assertEquals(ans[5].base_last, len(self.upstream) + len(self.front_bs) + 300 + 100)
        self.assertEquals(fragment_sequence[ans[5].base_first - 1:ans[5].base_last], "g" * 100)

        self.assertEquals(ans[6].feature.type, 'promoter')
        self.assertEquals(ans[6].feature.name, 'pBest')
        self.assertEquals(ans[6].feature.strand, -1)
        self.assertEquals(ans[6].base_first, len(self.upstream) + len(self.front_bs) + 400 + 1)
        self.assertEquals(ans[6].base_last, len(self.upstream) + len(self.front_bs) + 400 + 100)
        self.assertEquals(fragment_sequence[ans[6].base_first - 1:ans[6].base_last], "c" * 100)

    def test_passes_annotations_through_api(self):
        donor = "a" * 100 + "g" * 100 + "c" * 50 + "t" * 50 + "g" * 100 + "c" * 100
        cassette = ''.join([self.front_bs, donor, self.back_bs])
        flen = len(self.front_bs)

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)

        c = Genome.objects.get(pk=r['id'])
        f = c.fragments.all()[0].indexed_fragment()
        annotations = f.annotations()
        self.assertEquals(len(annotations), 7)

    def test_does_not_annotate_if_donor_has_overhangs(self):
        donor = "a" * 500
        cassette = ''.join(["(gg/)", self.front_bs, donor, self.back_bs])
        flen = len(self.front_bs)

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("does not have overhangs" in r["errors"], True)

    def test_does_not_annotate_if_annotation_is_missing_field(self):
        donor = "a" * 500
        cassette = ''.join([self.front_bs, donor, self.back_bs])
        flen = len(self.front_bs)

        annotations = [
          # missing base_first
          dict(base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("have all the required field" in r["errors"], True)

        annotations = [
          # missing base_last
          dict(base_first=flen + 1, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("have all the required field" in r["errors"], True)

        annotations = [
          # missing name
          dict(base_first=flen + 1, base_last=flen + 100, type="promoter",
               strand=1, qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("have all the required field" in r["errors"], True)

        annotations = [
          # missing type
          dict(base_first=flen + 1, base_last=flen + 100, name="pFavorite",
               strand=1, qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("have all the required field" in r["errors"], True)

        annotations = [
          # missing strand
          dict(base_first=flen + 1, base_last=flen + 100, name="pFavorite", type="promoter",
               qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
        r = json.loads(res.content)
        self.assertEquals("errors" in r, True)
        self.assertEquals("have all the required field" in r["errors"], True)

        annotations = [
          # not missing anything
          dict(base_first=flen + 1, base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
        ]

        res = self.client.post('/edge/genomes/%s/recombination/' % self.genome.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=self.arm_len,
                                                    create=True,
                                                    annotations=annotations,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)

    def test_adds_annotations_when_cassette_goes_in_reverse_direction(self):
        donor = "a" * 100 + "g" * 100 + "c" * 50 + "t" * 50 + "g" * 100 + "c" * 100
        cassette = ''.join([str(Seq(self.back_bs).reverse_complement()), donor,
                            str(Seq(self.front_bs).reverse_complement())])
        flen = len(self.back_bs)

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        c = recombine(self.genome, cassette, self.arm_len, annotations=annotations)
        f = c.fragments.all()[0].indexed_fragment()
        fragment_sequence = f.sequence

        ans = f.annotations()
        self.assertEquals(len(ans), 7)

        ans = sorted(ans, key=lambda a: (a.base_first, -a.base_last))

        self.assertEquals(ans[0].feature.type, 'operation')

        self.assertEquals(ans[1].feature.type, 'promoter')
        self.assertEquals(ans[1].feature.name, 'pBest')
        self.assertEquals(ans[1].feature.strand, 1)
        self.assertEquals(ans[1].base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(ans[1].base_last, len(self.upstream) + len(self.front_bs) + 100)
        self.assertEquals(fragment_sequence[ans[1].base_first - 1:ans[1].base_last], "g" * 100)

        self.assertEquals(ans[2].feature.type, 'gene')
        self.assertEquals(ans[2].feature.name, 'Best')
        self.assertEquals(ans[2].feature.strand, 1)
        self.assertEquals(ans[2].base_first, len(self.upstream) + len(self.front_bs) + 100 + 1)
        self.assertEquals(ans[2].base_last, len(self.upstream) + len(self.front_bs) + 100 + 100)
        self.assertEquals(fragment_sequence[ans[2].base_first - 1:ans[2].base_last], "c" * 100)

        self.assertEquals(ans[3].feature.type, 'terminator')
        self.assertEquals(ans[3].feature.name, 'tBest')
        self.assertEquals(ans[3].feature.strand, 1)
        self.assertEquals(ans[3].base_first, len(self.upstream) + len(self.front_bs) + 200 + 1)
        self.assertEquals(ans[3].base_last, len(self.upstream) + len(self.front_bs) + 200 + 50)
        self.assertEquals(fragment_sequence[ans[3].base_first - 1:ans[3].base_last], "a" * 50)

        self.assertEquals(ans[4].feature.type, 'terminator')
        self.assertEquals(ans[4].feature.name, 'tFavorite')
        self.assertEquals(ans[4].feature.strand, -1)
        self.assertEquals(ans[4].base_first, len(self.upstream) + len(self.front_bs) + 250 + 1)
        self.assertEquals(ans[4].base_last, len(self.upstream) + len(self.front_bs) + 250 + 50)
        self.assertEquals(fragment_sequence[ans[4].base_first - 1:ans[4].base_last], "g" * 50)

        self.assertEquals(ans[5].feature.type, 'gene')
        self.assertEquals(ans[5].feature.name, 'Favorite')
        self.assertEquals(ans[5].feature.strand, -1)
        self.assertEquals(ans[5].base_first, len(self.upstream) + len(self.front_bs) + 300 + 1)
        self.assertEquals(ans[5].base_last, len(self.upstream) + len(self.front_bs) + 300 + 100)
        self.assertEquals(fragment_sequence[ans[5].base_first - 1:ans[5].base_last], "c" * 100)

        self.assertEquals(ans[6].feature.type, 'promoter')
        self.assertEquals(ans[6].feature.name, 'pFavorite')
        self.assertEquals(ans[6].feature.strand, -1)
        self.assertEquals(ans[6].base_first, len(self.upstream) + len(self.front_bs) + 400 + 1)
        self.assertEquals(ans[6].base_last, len(self.upstream) + len(self.front_bs) + 400 + 100)
        self.assertEquals(fragment_sequence[ans[6].base_first - 1:ans[6].base_last], "t" * 100)

    def test_adds_annotations_correctly_when_homology_arm_has_extra_non_matching_bps(self):
        donor = "a" * 100 + "g" * 100 + "c" * 50 + "t" * 50 + "g" * 100 + "c" * 100
        cassette = ''.join(["aaa" + self.front_bs, donor, self.back_bs + "tt"])
        flen = len(self.front_bs) + 3

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        c = recombine(self.genome, cassette, self.arm_len, annotations=annotations)
        f = c.fragments.all()[0].indexed_fragment()
        fragment_sequence = f.sequence

        ans = f.annotations()
        self.assertEquals(len(ans), 7)

        ans = sorted(ans, key=lambda a: (a.base_first, -a.base_last))

        self.assertEquals(ans[0].feature.type, 'operation')

        self.assertEquals(ans[1].feature.type, 'promoter')
        self.assertEquals(ans[1].feature.name, 'pFavorite')
        self.assertEquals(ans[1].feature.strand, 1)
        self.assertEquals(ans[1].base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(ans[1].base_last, len(self.upstream) + len(self.front_bs) + 100)
        self.assertEquals(fragment_sequence[ans[1].base_first - 1:ans[1].base_last], "a" * 100)

        self.assertEquals(ans[2].feature.type, 'gene')
        self.assertEquals(ans[2].feature.name, 'Favorite')
        self.assertEquals(ans[2].feature.strand, 1)
        self.assertEquals(ans[2].base_first, len(self.upstream) + len(self.front_bs) + 100 + 1)
        self.assertEquals(ans[2].base_last, len(self.upstream) + len(self.front_bs) + 100 + 100)
        self.assertEquals(fragment_sequence[ans[2].base_first - 1:ans[2].base_last], "g" * 100)

        self.assertEquals(ans[3].feature.type, 'terminator')
        self.assertEquals(ans[3].feature.name, 'tFavorite')
        self.assertEquals(ans[3].feature.strand, 1)
        self.assertEquals(ans[3].base_first, len(self.upstream) + len(self.front_bs) + 200 + 1)
        self.assertEquals(ans[3].base_last, len(self.upstream) + len(self.front_bs) + 200 + 50)
        self.assertEquals(fragment_sequence[ans[3].base_first - 1:ans[3].base_last], "c" * 50)

        self.assertEquals(ans[4].feature.type, 'terminator')
        self.assertEquals(ans[4].feature.name, 'tBest')
        self.assertEquals(ans[4].feature.strand, -1)
        self.assertEquals(ans[4].base_first, len(self.upstream) + len(self.front_bs) + 250 + 1)
        self.assertEquals(ans[4].base_last, len(self.upstream) + len(self.front_bs) + 250 + 50)
        self.assertEquals(fragment_sequence[ans[4].base_first - 1:ans[4].base_last], "t" * 50)

        self.assertEquals(ans[5].feature.type, 'gene')
        self.assertEquals(ans[5].feature.name, 'Best')
        self.assertEquals(ans[5].feature.strand, -1)
        self.assertEquals(ans[5].base_first, len(self.upstream) + len(self.front_bs) + 300 + 1)
        self.assertEquals(ans[5].base_last, len(self.upstream) + len(self.front_bs) + 300 + 100)
        self.assertEquals(fragment_sequence[ans[5].base_first - 1:ans[5].base_last], "g" * 100)

        self.assertEquals(ans[6].feature.type, 'promoter')
        self.assertEquals(ans[6].feature.name, 'pBest')
        self.assertEquals(ans[6].feature.strand, -1)
        self.assertEquals(ans[6].base_first, len(self.upstream) + len(self.front_bs) + 400 + 1)
        self.assertEquals(ans[6].base_last, len(self.upstream) + len(self.front_bs) + 400 + 100)
        self.assertEquals(fragment_sequence[ans[6].base_first - 1:ans[6].base_last], "c" * 100)

    def test_adds_annotations_correctly_when_integrating_reversely_and_has_non_matching_bps(self):
        donor = "a" * 100 + "g" * 100 + "c" * 50 + "t" * 50 + "g" * 100 + "c" * 100
        cassette = ''.join(["aaa" + str(Seq(self.back_bs).reverse_complement()), donor,
                            str(Seq(self.front_bs).reverse_complement()) + "tt"])
        flen = len(self.back_bs) + 3

        annotations = [
          dict(base_first=flen + 1,   base_last=flen + 100, name="pFavorite", type="promoter",
               strand=1, qualifiers=None),
          dict(base_first=flen + 101, base_last=flen + 200, name="Favorite", type="gene",
               strand=1, qualifiers=dict(locus="Favorite", product="Favoritep")),
          dict(base_first=flen + 201, base_last=flen + 250, name="tFavorite", type="terminator",
               strand=1, qualifiers=None),
          dict(base_first=flen + 251, base_last=flen + 300, name="tBest", type="terminator",
               strand=-1, qualifiers=None),
          dict(base_first=flen + 301, base_last=flen + 400, name="Best", type="gene",
               strand=-1, qualifiers=dict(locus="Best", product="Bestp")),
          dict(base_first=flen + 401, base_last=flen + 500, name="pBest", type="promoter",
               strand=-1, qualifiers=None)
        ]

        c = recombine(self.genome, cassette, self.arm_len, annotations=annotations)
        f = c.fragments.all()[0].indexed_fragment()
        fragment_sequence = f.sequence

        ans = f.annotations()
        self.assertEquals(len(ans), 7)

        ans = sorted(ans, key=lambda a: (a.base_first, -a.base_last))

        self.assertEquals(ans[0].feature.type, 'operation')

        self.assertEquals(ans[1].feature.type, 'promoter')
        self.assertEquals(ans[1].feature.name, 'pBest')
        self.assertEquals(ans[1].feature.strand, 1)
        self.assertEquals(ans[1].base_first, len(self.upstream) + len(self.front_bs) + 1)
        self.assertEquals(ans[1].base_last, len(self.upstream) + len(self.front_bs) + 100)
        self.assertEquals(fragment_sequence[ans[1].base_first - 1:ans[1].base_last], "g" * 100)

        self.assertEquals(ans[2].feature.type, 'gene')
        self.assertEquals(ans[2].feature.name, 'Best')
        self.assertEquals(ans[2].feature.strand, 1)
        self.assertEquals(ans[2].base_first, len(self.upstream) + len(self.front_bs) + 100 + 1)
        self.assertEquals(ans[2].base_last, len(self.upstream) + len(self.front_bs) + 100 + 100)
        self.assertEquals(fragment_sequence[ans[2].base_first - 1:ans[2].base_last], "c" * 100)

        self.assertEquals(ans[3].feature.type, 'terminator')
        self.assertEquals(ans[3].feature.name, 'tBest')
        self.assertEquals(ans[3].feature.strand, 1)
        self.assertEquals(ans[3].base_first, len(self.upstream) + len(self.front_bs) + 200 + 1)
        self.assertEquals(ans[3].base_last, len(self.upstream) + len(self.front_bs) + 200 + 50)
        self.assertEquals(fragment_sequence[ans[3].base_first - 1:ans[3].base_last], "a" * 50)

        self.assertEquals(ans[4].feature.type, 'terminator')
        self.assertEquals(ans[4].feature.name, 'tFavorite')
        self.assertEquals(ans[4].feature.strand, -1)
        self.assertEquals(ans[4].base_first, len(self.upstream) + len(self.front_bs) + 250 + 1)
        self.assertEquals(ans[4].base_last, len(self.upstream) + len(self.front_bs) + 250 + 50)
        self.assertEquals(fragment_sequence[ans[4].base_first - 1:ans[4].base_last], "g" * 50)

        self.assertEquals(ans[5].feature.type, 'gene')
        self.assertEquals(ans[5].feature.name, 'Favorite')
        self.assertEquals(ans[5].feature.strand, -1)
        self.assertEquals(ans[5].base_first, len(self.upstream) + len(self.front_bs) + 300 + 1)
        self.assertEquals(ans[5].base_last, len(self.upstream) + len(self.front_bs) + 300 + 100)
        self.assertEquals(fragment_sequence[ans[5].base_first - 1:ans[5].base_last], "c" * 100)

        self.assertEquals(ans[6].feature.type, 'promoter')
        self.assertEquals(ans[6].feature.name, 'pFavorite')
        self.assertEquals(ans[6].feature.strand, -1)
        self.assertEquals(ans[6].base_first, len(self.upstream) + len(self.front_bs) + 400 + 1)
        self.assertEquals(ans[6].base_last, len(self.upstream) + len(self.front_bs) + 400 + 100)
        self.assertEquals(fragment_sequence[ans[6].base_first - 1:ans[6].base_last], "t" * 100)
