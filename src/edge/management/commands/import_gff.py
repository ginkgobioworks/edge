from django.db import transaction
from django.core.management.base import BaseCommand
from edge.models import Genome

class Command(BaseCommand):

    @transaction.atomic
    def handle(self, *args, **options):
        if len(args) != 2:
            raise Exception('Expecting two arguments: name of genome and GFF file')
        g = Genome.create(args[0])
        u = g.edit()
        u.import_gff(args[1])
