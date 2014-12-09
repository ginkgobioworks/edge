import os
import json
from django.test import TestCase
from edge.models import Genome, Operation, Fragment, Genome_Fragment
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn
from edge.crispr import find_crispr_target
from Bio.Seq import Seq


class GenomeCrisprDSBTest(TestCase):

    def build_genome(self, *sequences):
        g = Genome(name='Foo')
        g.save()
        for seq in sequences:
            f = Fragment.create_with_sequence('Bar', seq, circular=True)
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
        g = self.build_genome(s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].target_start, s1.index(guide)+1)
        self.assertEquals(t[0].target_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_target_on_reverse_strand(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome(str(Seq(s1+pam+s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].target_end, len(s2)+3+1)
        self.assertEquals(t[0].target_start, len(s2)+3+1+20-1)
        self.assertEquals(t[0].pam, 'ngg')

    def test_find_crispr_target_finds_multiple_crispr_targets(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        g = self.build_genome((s1+pam+s2)+str(Seq(s1+pam+s2).reverse_complement()))
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 2)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].target_start, s1.index(guide)+1)
        self.assertEquals(t[0].target_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')
        self.assertEquals(t[1].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[1].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[1].target_end, len(s1+pam+s2)+len(s2)+3+1)
        self.assertEquals(t[1].target_start, len(s1+pam+s2)+len(s2)+3+1+20-1)
        self.assertEquals(t[1].pam, 'ngg')

    def test_find_crispr_target_does_not_find_target_without_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'ccc'
        g = self.build_genome(s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_does_not_find_target_with_part_of_pam(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgc'
        g = self.build_genome(s1+pam+s2)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 0)

    def test_find_crispr_target_finds_target_across_circular_boundary(self):
        s1 = 'agaaggtctggtagcgatgtagtcgatct'
        s2 = 'gactaggtacgtagtcgtcaggtcagtca'
        pam = 'cgg'
        s = s1[-10:]+pam+s2+s1[0:-10]
        print s
        print s1+pam+s2
        g = self.build_genome(s)
        guide = s1[-20:]
        t = find_crispr_target(g, guide, 'ngg')
        self.assertEquals(len(t), 1)
        self.assertEquals(t[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(t[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(t[0].target_start, s1.index(guide)+1)
        self.assertEquals(t[0].target_end, len(s1))
        self.assertEquals(t[0].pam, 'ngg')

    def test_crispr_dsb_creates_operations(self):
        self.assertEquals(True, False)

    def test_crispr_dsb_finds_and_annotates_target_on_forward_strand(self):
        self.assertEquals(True, False)

    def test_crispr_dsb_finds_and_annotates_target_on_reverse_strand(self):
        self.assertEquals(True, False)

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
