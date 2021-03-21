from celery import shared_task
from edge.models import Genome, Operation
from edge.blastdb import build_genome_db
from edge.recombine import annotate_integration
from django.db import transaction


@shared_task
def build_genome_blastdb(genome_id):
    genome = Genome.objects.get(pk=genome_id)
    build_genome_db(genome, refresh=True)


@transaction.atomic()
def annotate_integration_on_genome(genome, new_genome,
                                   regions_before, regions_after,
                                   cassette_name, op, cassette, annotations):
    annotate_integration(genome, new_genome, regions_before, regions_after,
                         cassette_name, op, cassette, annotations)


@shared_task
def annotate_integration_task(genome_id, new_genome_id,
                              regions_before, regions_after,
                              cassette_name, op_id, cassette, annotations):
    genome = Genome.objects.get(pk=genome_id)
    new_genome = Genome.objects.get(pk=new_genome_id)
    op = Operation.objects.get(pk=op_id)
    annotate_integration_on_genome(genome, new_genome,
                                   regions_before, regions_after,
                                   cassette_name, op, cassette, annotations)
