import os
import json
from Bio.Seq import Seq
from django.test import TestCase


class GenomeBlastTest(TestCase):

    def test_blast_finds_sequence_on_specified_genome(self):
        from edge.models import Genome, Fragment, Genome_Fragment
        from edge.management.commands.build_edge_blastdb import build_fragment_db, fragment_fasta_fn

        s1 = 'atcggtatcttctatgcgtatgcgtcatgattatatatattagcggcatg'
        s2 = 'agcgtcgatgcatgagtcgatcggcagtcgtgtagtcgtcgtatgcgtta'
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
        except:
            pass
        build_fragment_db()

        query = s1[6:20]+'aaaaaaaaa'

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
        from edge.models import Genome, Fragment, Genome_Fragment
        from edge.management.commands.build_edge_blastdb import build_fragment_db, fragment_fasta_fn

        s1 = 'atcggtatcttctatgcgtatgcgtcatgattatatatattagcggcatg'
        g1 = Genome(name='Foo')
        g1.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        try:
            os.unlink(fragment_fasta_fn(f1))
        except:
            pass
        build_fragment_db()

        query = str(Seq(s1[6:20]).reverse_complement())+'tttttttttt'

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
        self.assertEquals(d[0]['subject_start'], 20)
        self.assertEquals(d[0]['subject_end'], 7)
