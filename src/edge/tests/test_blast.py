import os
import json
import random

from Bio.Seq import Seq
from django.test import TestCase

from edge import get_random_sequence
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn
from edge.blast import blast_genome
from edge.models import (
    Fragment,
    Genome,
    Genome_Fragment,
)


class GenomeBlastTest(TestCase):
    def setUp(self):
        random.seed(0)

    def test_finds_sequence_on_specified_genome(self):
        s1 = get_random_sequence(200)
        s2 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        f2 = Fragment.create_with_sequence('Baz', s2)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()
        Genome_Fragment(genome=g1, fragment=f2, inherited=False).save()

        g2 = Genome(name='Far')
        g2.save()
        f3 = Fragment.create_with_sequence('Bar', s1)
        Genome_Fragment(genome=g2, fragment=f3, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
            os.unlink(fragment_fasta_fn(f2))
            os.unlink(fragment_fasta_fn(f3))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)
        g1 = Genome.objects.get(pk=g1.id)

        query = s1[6:20] + 'aaaaaaaaa'
        r = blast_genome(g1, 'blastn', query)
        # only returns hit from genome
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, f1.id)
        self.assertEquals(r[0].query_start, 1)
        self.assertEquals(r[0].query_end, 14)
        self.assertEquals(r[0].subject_start, 7)
        self.assertEquals(r[0].subject_end, 20)
        self.assertEquals(r[0].strand(), 1)

    def test_aligns_sequence_to_antisense_strand(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)
        g1 = Genome.objects.get(pk=g1.id)

        query = str(Seq(s1[6:20]).reverse_complement()) + 'tttttttttt'
        r = blast_genome(g1, 'blastn', query)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, f1.id)
        self.assertEquals(r[0].query_start, 1)
        self.assertEquals(r[0].query_end, 14)
        self.assertEquals(r[0].subject_start, 20)
        self.assertEquals(r[0].subject_end, 7)
        self.assertEquals(r[0].strand(), -1)

    def test_aligns_sequence_across_boundry_for_circular_fragment(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1, circular=True)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)
        g1 = Genome.objects.get(pk=g1.id)

        query = (s1[-10:] + s1[0:10]) + 'ttttttttttt'
        res = blast_genome(g1, 'blastn', query)

        # we are not removing redundant matches when matching across circular
        # boundaries, since blasting across circular boundary of a genome is a
        # rare case. so in this particular case, you will find two results, one
        # for the end of the query at the start of the genome, one for across
        # the circular boundary.

        found = False
        for r in res:
            if r.query_start == 1 and r.query_end == 20:
                self.assertEquals(r.fragment_id, f1.id)
                self.assertEquals(r.query_start, 1)
                self.assertEquals(r.query_end, 20)
                self.assertEquals(r.subject_start, len(s1) - 10 + 1)
                self.assertEquals(r.subject_end, len(s1) + 10)
                self.assertEquals(r.fragment_length, len(s1))
                self.assertEquals(r.strand(), 1)
                found = True
                break
        self.assertEquals(found, True)

    def test_does_not_return_duplicate_hits_for_circular_fragments(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1, circular=True)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)
        g1 = Genome.objects.get(pk=g1.id)

        query = s1[5:20] + 'tttttttttt'
        r = blast_genome(g1, 'blastn', query)
        self.assertEquals(len(r), 1)

    def test_does_not_align_sequence_across_boundry_for_non_circular_fragment(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1, circular=False)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)
        g1 = Genome.objects.get(pk=g1.id)

        query = (s1[-10:] + s1[0:10]) + 'tttttttttt'
        res = blast_genome(g1, 'blastn', query)

        for r in res:
            self.assertEquals(r.subject_start > 0 and r.subject_start <= len(s1), True)
            self.assertEquals(r.subject_end > 0 and r.subject_end <= len(s1), True)


class GenomeBlastAPITest(TestCase):
    def setUp(self):
        random.seed(0)

    def test_blast_finds_sequence_on_specified_genome(self):
        s1 = get_random_sequence(200)
        s2 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        f2 = Fragment.create_with_sequence('Baz', s2)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()
        Genome_Fragment(genome=g1, fragment=f2, inherited=False).save()

        g2 = Genome(name='Far')
        g2.save()
        f3 = Fragment.create_with_sequence('Bar', s1)
        Genome_Fragment(genome=g2, fragment=f3, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
            os.unlink(fragment_fasta_fn(f2))
            os.unlink(fragment_fasta_fn(f3))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)

        query = s1[6:20] + 'aaaaaaaaa'

        res = self.client.post('/edge/genomes/%s/blast/' % g1.id,
                               data=json.dumps(dict(program='blastn', query=query)),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)

        # only returns hit from genome
        self.assertEquals(len(d), 1)
        self.assertEquals(d[0]['fragment_id'], f1.id)
        self.assertEquals(d[0]['query_start'], 1)
        self.assertEquals(d[0]['query_end'], 14)
        self.assertEquals(d[0]['subject_start'], 7)
        self.assertEquals(d[0]['subject_end'], 20)

        # blast in other genome works too
        res = self.client.post('/edge/genomes/%s/blast/' % g2.id,
                               data=json.dumps(dict(program='blastn', query=query)),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(len(d), 1)
        self.assertEquals(d[0]['fragment_id'], f3.id)

    def test_blast_aligns_sequence_to_antisense_strand(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except BaseException:
            pass
        build_all_genome_dbs(refresh=True)

        query = str(Seq(s1[6:20]).reverse_complement()) + 'tttttttttt'

        res = self.client.post('/edge/genomes/%s/blast/' % g1.id,
                               data=json.dumps(dict(program='blastn', query=query)),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)

        self.assertEquals(len(d), 1)
        self.assertEquals(d[0]['fragment_id'], f1.id)
        self.assertEquals(d[0]['query_start'], 1)
        self.assertEquals(d[0]['query_end'], 14)
        self.assertEquals(d[0]['subject_start'], 20)
        self.assertEquals(d[0]['subject_end'], 7)
