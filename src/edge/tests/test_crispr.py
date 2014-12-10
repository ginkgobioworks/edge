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
            except:
                pass
        build_all_genome_dbs(refresh=True)
        return g

    def test_find_crispr_target_finds_target_on_forward_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, s1.index(guide)+1)
        self.assertEquals(t[0].subject_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_only_finds_perfect_match_to_guide(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1+pam+s2)
        guide = 'aaaaa'+s1[-15:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_finds_target_on_reverse_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, str(Seq(s1+pam+s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_end, len(s2)+3+1)
        self.assertEquals(t[0].subject_start, len(s2)+3+1+20-1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_multiple_crispr_targets(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, (s1+pam+s2)+str(Seq(s1+pam+s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 2)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, s1.index(guide)+1)
        self.assertEquals(t[0].subject_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')
        self.assertEquals(t[1].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[1].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[1].subject_end, len(s1+pam+s2)+len(s2)+3+1)
        self.assertEquals(t[1].subject_start, len(s1+pam+s2)+len(s2)+3+1+20-1)
        self.assertEquals(t[1].pam, 'ngg')

    def test_find_crispr_target_does_not_find_target_without_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'ccc'
        g = self.build_genome(False, s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_does_not_find_target_with_part_of_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgc'
        g = self.build_genome(False, s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_finds_target_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[10:]+pam+s2+s1[0:10]
        g = self.build_genome(True, s)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, (s1.index(guide)+1-10-1)%len(s)+1)
        self.assertEquals(t[0].subject_end, (len(s1)-10-1)%len(s)+1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_target_with_pam_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = pam[1:]+s2+s1+pam[:1]
        g = self.build_genome(True, s)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, len(pam[1:]+s2)+s1.index(guide)+1)
        self.assertEquals(t[0].subject_end, len(s)-1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_reverse_complement_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[10:]+pam+s2+s1[0:10]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, (len(s)-(s1.index(guide)+1-10))%len(s)+1)
        self.assertEquals(t[0].subject_end, (len(s)-(len(s1)-10))%len(s)+1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_reverse_complement_with_pam_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = pam[1:]+s2+s1+pam[:1]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_end, 2)
        self.assertEquals(t[0].subject_start, 2+len(guide)-1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_crispr_dsb_finds_and_annotates_target_on_forward_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1+pam+s2)
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg', 3)
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = g.fragments.all()[0].indexed_fragment().annotations()
        # XXX
        # self.assertEquals(len(a), 0)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, s1.index(guide)+1)
        self.assertEquals(a[0].base_last, len(s1))
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(guide))
        self.assertEquals(a[0].feature.strand, 1)
        self.assertEquals(a[0].feature.name, 'CRISPR DSB target')

    def test_crispr_dsb_finds_and_annotates_target_on_reverse_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, str(Seq(s1+pam+s2).reverse_complement()))
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg', 3)
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = g.fragments.all()[0].indexed_fragment().annotations()
        # XXX
        # self.assertEquals(len(a), 0)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, len(s2)+3+1)
        self.assertEquals(a[0].base_last, len(s2)+3+1+20-1)
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(guide))
        self.assertEquals(a[0].feature.strand, -1)
        self.assertEquals(a[0].feature.name, 'CRISPR DSB target')

    def test_crispr_dsb_finds_and_annotates_target_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[-15:]+pam+s2+s1[0:len(s1)-15]
        g = self.build_genome(True, s)
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg', 3)
        self.assertNotEquals(c.id, g.id)
        self.assertEquals(c.parent.id, g.id)

        a = g.fragments.all()[0].indexed_fragment().annotations()
        # XXX
        # self.assertEquals(len(a), 0)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 2)
        self.assertEquals(a[1].base_first, len(s)-5+1)
        self.assertEquals(a[1].base_last, len(s))
        self.assertEquals(a[1].feature_base_first, 1)
        self.assertEquals(a[1].feature_base_last, 5)
        self.assertEquals(a[1].feature.strand, 1)
        self.assertEquals(a[1].feature.name, 'CRISPR DSB target')
        self.assertEquals(True, False)

    def test_crispr_dsb_finds_and_annotates_reverse_complement_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[10:]+pam+s2+s1[0:10]
        g = self.build_genome(True, str(Seq(s).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].subject_start, (len(s)-(s1.index(guide)+1-10))%len(s)+1)
        self.assertEquals(t[0].subject_end, (len(s)-(len(s1)-10))%len(s)+1)
        self.assertEquals(t[0].pam, 'ngg')
        self.assertEquals(True, False)

    def test_crispr_dsb_creates_operations(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(False, s1+pam+s2)
        guide = s1[-20:]

        c = crispr_dsb(g, guide, 'ngg', 3)
        self.assertEquals(g.operations.count(), 0)
        self.assertEquals(c.operations.count(), 1)
        self.assertEquals(c.operations.all()[0].type, Operation.CRISPR_DSB[0])

    def test_crispr_dsb_api_works(self):
        self.assertEquals(Operation.objects.count(), 0)
        data = dict(genome_name='FooBar', notes='blah',
                    guide=self.target, pam='ngg', create=True)
        res = self.client.post('/genomes/'+str(self.genome.id)+'/crispr/dsb/',
                               data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(Operation.objects.count(), 1)
        self.assertEquals(self.genome.children.count(), 1)
        c = self.genome.children.all()[0]
        self.assertEquals(c.operations.all()[0].type, Operation.CRISPR_DSB[0])
        self.assertEquals(c.params, json.dumps(dict(guide=self.target, pam='ngg')))
        self.assertEquals(True, False)
