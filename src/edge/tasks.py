from celery import task
from edge.models import Genome
from edge.blastdb import build_genome_db

@task(name='build_genome_blastdb')
def build_genome_blastdb(genome_id):
  genome = Genome.objects.get(pk=genome_id)
  build_genome_db(genome, refresh=True)
