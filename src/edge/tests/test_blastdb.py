import os
import random
from Bio import SeqIO
from django.test import TestCase
from edge import get_random_sequence
from edge.blastdb import build_genome_db, build_fragment_fasta, fragment_fasta_fn
from edge.models import (
    Fragment,
    Genome,
    Genome_Fragment,
)


class BuildBlastDBTest(TestCase):
    def setUp(self):
        random.seed(0)

    def test_builds_fragment_fastas(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name="Foo")
        g1.save()
        f1 = Fragment.create_with_sequence("Bar", s1)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        fn = fragment_fasta_fn(f1)
        try:
            os.unlink(fn)
        except BaseException:
            pass

        fn = build_fragment_fasta(f1)
        records = list(SeqIO.parse(fn, "fasta"))
        self.assertEquals(len(records), 1)
        self.assertEquals(str(records[0].seq), s1)

    def test_builds_genome_db_with_different_names_on_separate_attempts(self):
        s1 = get_random_sequence(200)
        g1 = Genome(name="Foo")
        g1.save()
        f1 = Fragment.create_with_sequence("Bar", s1)
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()

        dbname1 = build_genome_db(g1)
        self.assertEquals(dbname1 is None, False)
        self.assertEquals(dbname1, g1.blastdb)

        g1.blastdb = None
        dbname2 = build_genome_db(g1)
        self.assertEquals(dbname2 == dbname1, False)
