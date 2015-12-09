from django.core.management.base import BaseCommand
from django.db import transaction
from edge.models import Genome
from edge.management.commands.remove_fragment import remove_fragment


@transaction.atomic()
def remove_genome(genome_id):
    genome = Genome.objects.get(pk=genome_id)
    fragments = [f.id for f in genome.fragments.all()]
    for fid in fragments:
        remove_fragment(fid)
    genome.delete()


class Command(BaseCommand):

    def handle(self, *args, **options):
        if len(args) != 1:
            raise Exception('Expecting genome ID as argument')
        remove_genome(args[0])
