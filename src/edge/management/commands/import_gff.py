from django.core.management.base import BaseCommand
from edge import import_gff


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('genome', type=str)
        parser.add_argument('filename', type=str)

    def handle(self, *args, **options):
        import_gff(options['genome'], options['filename'])
