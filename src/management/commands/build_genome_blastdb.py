from edge.blastdb import build_all_genome_dbs
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    def handle(self, *args, **options):
        build_all_genome_dbs()
