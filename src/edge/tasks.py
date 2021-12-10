from celery import shared_task
from edge.models import Genome
from edge.blastdb import build_genome_db


@shared_task
def build_genome_fragment_indices(genome_id):
    genome = Genome.objects.get(pk=genome_id)
    genome.indexed_genome()


@shared_task
def build_genome_blastdb(genome_id):
    genome = Genome.objects.get(pk=genome_id)
    build_genome_db(genome, refresh=True)
