from django.core.management.base import BaseCommand
from django.db import transaction
from edge.models import Fragment, Chunk_Feature, Chunk


@transaction.atomic()
def remove_fragment(fragment_id):
    fragment = Fragment.objects.get(pk=fragment_id)
    fragment.fragment_chunk_location_set.all().delete()
    fragment.edge_set.all().delete()
    Chunk_Feature.objects.filter(chunk__initial_fragment_id=fragment_id).delete()
    fragment.start_chunk = None
    fragment.save()
    Chunk.objects.filter(initial_fragment_id=fragment_id).delete()
    fragment.genome_fragment_set.all().delete()
    fragment.delete()


class Command(BaseCommand):

    def handle(self, *args, **options):
        if len(args) != 1:
            raise Exception('Expecting fragment ID as argument')
        remove_fragment(args[0])
