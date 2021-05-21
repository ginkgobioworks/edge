# flake8: noqa
# imports multiple GFF+FASTA files
#
# python import_multi.py genome_name <file names...>
#
#
# Example:
# python split_gff.py GRCh38 GRCh38_latest_genomic.gff GRCh38_latest_genomic.fna
# in docker, run
# time scripts/_python example/import_multi.py 'GRCh38.2' example/GRCh38/*combined.gff
# time scripts/_python example/import_multi.py GCA_003668045.2 example/GCA_003668045.2/*combined.gff

import django

django.setup()

import sys
from edge.models import Genome
from edge.importer import GFFImporter

name = sys.argv[1]
fns = sys.argv[2:]

"""
if Genome.objects.filter(name=name).count() > 0:
  raise Exception('There is already a genome named "%s"' % (name,))
"""

if Genome.objects.filter(name=name).count() > 0:
    genome = Genome.objects.get(name=name)
else:
    genome = Genome.create(name)

for gff_fasta_fn in fns:
    print(gff_fasta_fn)
    GFFImporter(genome, gff_fasta_fn).do_import()
