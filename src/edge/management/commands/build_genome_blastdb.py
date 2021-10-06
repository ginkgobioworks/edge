from edge.blastdb import build_all_genome_dbs
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            '--refresh',
            action='store_true',
            help='Rebuild BLAST database files',
        )

    def handle(self, *args, **options):
        if options['refresh']:
            build_all_genome_dbs(refresh=True)
        else:
            build_all_genome_dbs(refresh=False)
