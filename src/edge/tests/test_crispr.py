import os
import json
from django.test import TestCase
from edge.models import Genome, Operation, Fragment, Genome_Fragment
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn
from edge.crispr import find_crispr_target, crispr_dsb
from Bio.Seq import Seq


class GenomeCrisprDSBTest(TestCase):

    def build_genome(self, circular, *sequences):
        g = Genome(name='Foo')
        g.save()
        for seq in sequences:
            f = Fragment.create_with_sequence('Bar', seq, circular=circular)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
            try:
                os.unlink(fragment_fasta_fn(f))
            except OSError:
                pass
        build_all_genome_dbs(refresh=True)
        return Genome.objects.get(pk=g.id)

    def test_find_crispr_target_finds_target_on_forward_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, s1.index(guide) + 1)
        self.assertEquals(t[0].subject_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_only_finds_perfect_match_to_guide(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = 'aaaaa' + s1[-15:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_finds_target_on_reverse_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, str(Seq(s1 + pam + s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_end, len(s2) + 3 + 1)
        self.assertEquals(t[0].subject_start, len(s2) + 3 + 1 + 20 - 1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_multiple_crispr_targets(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, (s1 + pam + s2) + str(Seq(s1 + pam + s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 2)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, s1.index(guide) + 1)
        self.assertEquals(t[0].subject_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')
        self.assertEquals(t[1].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[1].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[1].subject_end, len(s1 + pam + s2) + len(s2) + 3 + 1)
        self.assertEquals(t[1].subject_start, len(s1 + pam + s2) + len(s2) + 3 + 1 + 20 - 1)
        self.assertEquals(t[1].pam, 'ngg')

    def test_find_crispr_target_does_not_find_target_without_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'ccc'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_does_not_find_target_with_part_of_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgc'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_finds_target_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[10:] + pam + s2 + s1[0:10]
        g = self.build_genome(True, s)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, (s1.index(guide) + 1 - 10 - 1) % len(s) + 1)
        self.assertEquals(t[0].subject_end, (len(s1) - 10 - 1) % len(s) + 1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_target_with_pam_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = pam[1:] + s2 + s1 + pam[:1]
        g = self.build_genome(True, s)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, len(pam[1:] + s2) + s1.index(guide) + 1)
        self.assertEquals(t[0].subject_end, len(s) - 1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_reverse_complement_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[10:] + pam + s2 + s1[0:10]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, (len(s) - (s1.index(guide) + 1 - 10)) % len(s) + 1)
        self.assertEquals(t[0].subject_end, (len(s) - (len(s1) - 10)) % len(s) + 1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_reverse_complement_with_pam_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = pam[1:] + s2 + s1 + pam[:1]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_end, 2)
        self.assertEquals(t[0].subject_start, 2 + len(guide) - 1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_crispr_dsb_finds_and_annotates_target_on_forward_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]

        a = g.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 0)

        c = crispr_dsb(g, guide, 'ngg')
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, s1.index(guide) + 1)
        self.assertEquals(a[0].base_last, len(s1))
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(guide))
        self.assertEquals(a[0].feature.strand, 1)
        self.assertEquals(a[0].feature.name, 'CRISPR-Cas9 (pam ngg) target')
        self.assertEquals(a[0].feature.operation.type, Operation.CRISPR_DSB[0])
        self.assertEquals(a[0].feature.operation.genome, c)

        # annotation is visible in parent genome, since it's not on new base pairs
        a = g.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].feature.operation.genome, c)

    def test_crispr_dsb_finds_and_annotates_target_on_reverse_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, str(Seq(s1 + pam + s2).reverse_complement()))
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg')
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, len(s2) + 3 + 1)
        self.assertEquals(a[0].base_last, len(s2) + 3 + 1 + 20 - 1)
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(guide))
        self.assertEquals(a[0].feature.strand, -1)
        self.assertEquals(a[0].feature.name, 'CRISPR-Cas9 (pam ngg) target')
        self.assertEquals(a[0].feature.operation.type, Operation.CRISPR_DSB[0])
        self.assertEquals(a[0].feature.operation.genome, c)

    def test_crispr_dsb_finds_and_annotates_target_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[-15:] + pam + s2 + s1[0:len(s1) - 15]
        g = self.build_genome(True, s)
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg')
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 2)
        self.assertEquals(a[0].base_first, 1)
        self.assertEquals(a[0].base_last, 15)
        self.assertEquals(a[0].feature_base_first, 6)
        self.assertEquals(a[0].feature_base_last, 20)
        self.assertEquals(a[0].feature.strand, 1)
        self.assertEquals(a[1].base_first, len(s) - 5 + 1)
        self.assertEquals(a[1].base_last, len(s))
        self.assertEquals(a[1].feature_base_first, 1)
        self.assertEquals(a[1].feature_base_last, 5)
        self.assertEquals(a[1].feature.strand, 1)

    def test_crispr_dsb_finds_and_annotates_reverse_complement_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[-15:] + pam + s2 + s1[0:len(s1) - 15]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg')
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 2)
        self.assertEquals(a[0].base_first, 1)
        self.assertEquals(a[0].base_last, 5)
        self.assertEquals(a[0].feature_base_first, 16)
        self.assertEquals(a[0].feature_base_last, 20)
        self.assertEquals(a[0].feature.strand, -1)

        self.assertEquals(a[1].base_first, len(s) - 15 + 1)
        self.assertEquals(a[1].base_last, len(s))
        self.assertEquals(a[1].feature_base_first, 1)
        self.assertEquals(a[1].feature_base_last, 15)
        self.assertEquals(a[1].feature.strand, -1)

    def test_crispr_dsb_creates_operations(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg')
        self.assertEquals(g.operation_set.count(), 0)
        self.assertEquals(c.operation_set.count(), 1)
        self.assertEquals(c.operation_set.all()[0].type, Operation.CRISPR_DSB[0])

    def test_crispr_dsb_creates_new_fragment(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]
        c = crispr_dsb(g, guide, 'ngg')
        self.assertEquals(c.fragments.all()[0].parent, g.fragments.all()[0])

    def test_crispr_dsb_api_works(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]

        data = dict(genome_name='FooBar', notes='blah', guide=guide, pam='ngg', create=True)
        res = self.client.post('/edge/genomes/' + str(g.id) + '/crispr/dsb/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)

        self.assertEquals(Operation.objects.count(), 1)
        self.assertEquals(g.children.count(), 1)
        c = g.children.all()[0]
        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].feature.operation.genome, c)
        self.assertEquals(c.operation_set.all()[0].type, Operation.CRISPR_DSB[0])
        self.assertEquals(c.operation_set.all()[0].params, json.dumps(dict(guide=guide, pam='ngg')))

    def test_multiple_api_calls_return_same_child(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1 + pam + s2)
        guide = s1[-20:]
        data = dict(genome_name='FooBar', notes='blah', guide=guide, pam='ngg', create=True)

        res = self.client.post('/edge/genomes/' + str(g.id) + '/crispr/dsb/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)
        c1 = r['id']

        res = self.client.post('/edge/genomes/' + str(g.id) + '/crispr/dsb/', data=json.dumps(data),
                               content_type='application/json')
        # returns 200 not 201
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)
        c2 = r['id']
        self.assertEquals(c1, c2)
