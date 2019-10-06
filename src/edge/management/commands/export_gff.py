from django.core.management.base import BaseCommand
from edge.models.genome import Genome
from edge.io import IO


class Command(BaseCommand):

    def handle(self, *args, **options):
        if len(args) != 2:
            raise Exception('Expecting integer genome ID and filename as arguments')
        try:
            genome_id = int(args[0])
        except BaseException:
            raise Exception('Expecting integer genome ID and filename as arguments')

        io = IO(Genome.objects.get(pk=genome_id))
        io.to_gff(args[1])
