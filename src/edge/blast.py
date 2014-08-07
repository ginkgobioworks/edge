from django.conf import settings
from edge.models import Fragment


BLAST_DB = "%s/edge-nucl" % settings.NCBI_DATA_DIR


class Blast_Accession(object):

    @staticmethod
    def make(fragment):
        return '%s' % fragment.id

    def __init__(self, accession):
        self.fragment_id = int(accession)

    @property
    def fragment(self):
        return Fragment.objects.get(pk=self.fragment_id)
