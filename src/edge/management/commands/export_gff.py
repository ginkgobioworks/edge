from django.core.management.base import BaseCommand
from edge.models.genome import Genome
from edge.io import IO


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('genome-id', type=int)
        parser.add_argument('filename', type=str)

    def handle(self, *args, **options):
        io = IO(Genome.objects.get(pk=options['genome-id']))
        io.to_gff(options['filename'])
