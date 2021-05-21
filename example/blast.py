# flake8: noqa
import django

django.setup()

import sys
from edge.models import Genome
from edge.blast import blast_genome
from edge.pcr import pcr_from_genome

genome = Genome.objects.get(pk=36161)
print(genome.fragments.count())

res = blast_genome(genome, "blastn", "ctcacgttactgtggggtggaggggaca")
for r in res:
    print("p: %s" % r.evalue)

primer_fwd = "aggaagtgccattccgcctgacctcgtctcactgaccgtctctctcctgagtccgga"
primer_rev = "aaagtgtcaaggtctcacgttactgtggggtggaggggaca"

res = blast_genome(genome, "blastn", primer_fwd)
for r in res:
    print("fwd: %s" % r.evalue)

res = blast_genome(genome, "blastn", primer_rev)
for r in res:
    print("rev: %s" % r.evalue)

print(pcr_from_genome(genome, primer_fwd, primer_rev))
