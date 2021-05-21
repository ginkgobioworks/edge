# flake8: noqa
import django

django.setup()

import sys
from edge.models import Genome
from edge.blastdb import build_genome_db

genome = Genome.objects.get(name="GRCh38.3")
assert genome.fragments.count() == 639

build_genome_db(genome)
