from django.core.management.base import BaseCommand
from edge import import_gff


class Command(BaseCommand):

    def handle(self, *args, **options):
        if len(args) != 2:
            raise Exception('Expecting two arguments: name of genome and GFF file')
        import_gff(args[0], args[1])
