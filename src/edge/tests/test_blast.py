import os
import json
from django.test import TestCase


class GenomeBlastTest(TestCase):

    def test_blast_finds_sequence(self):
        from edge.models import Genome, Fragment, Genome_Fragment
        from edge.management.commands.build_blastdb import build_fragment_db, fragment_fasta_fn

        s1 = 'atcggtatcttctatgcgtatgcgtcatgattatatatattagcggcatg'
        s2 = 'agcgtcgatgcatgagtcgatcggcagtcgtgtagtcgtcgtatgcgtta'
        g = Genome(name='Foo')
        g.save()
        f1 = Fragment.create_with_sequence('Bar', s1)
        f2 = Fragment.create_with_sequence('Baz', s2)
        Genome_Fragment(genome=g, fragment=f1, inherited=False).save()
        Genome_Fragment(genome=g, fragment=f2, inherited=False).save()

        os.unlink(fragment_fasta_fn(f1))
        os.unlink(fragment_fasta_fn(f2))
        build_fragment_db()

        res = self.client.post('/edge/genomes/%s/blast/' % g.id,
                               data=json.dumps(dict(program='blastn', query=s1[6:20])),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(len(d), 1)
        self.assertEquals(d[0]['fragment_id'], f1.id)
        self.assertEquals(d[0]['query_start'], 1)
        self.assertEquals(d[0]['query_end'], 14)
        self.assertEquals(d[0]['subject_start'], 7)
        self.assertEquals(d[0]['subject_end'], 20)
