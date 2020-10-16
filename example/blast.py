import django
django.setup()

import sys
from edge.models import Genome
from edge.blast import blast_genome
from edge.pcr import pcr_from_genome

genome = Genome.objects.get(pk=36161)
print(genome.fragments.count())

primer = "cacgttactgtggggtggaggggacag"
print(blast_genome(genome, "blastn", primer))

primer_fwd = "aggaagtgccattccgcctgacctcgtctcactgaccgtctctctcctgagtccgga"
primer_rev = "aaagtgtcaaggtctcacgttactgtggggtggaggggacag"
print(pcr_from_genome(genome, primer_fwd, primer_rev))
