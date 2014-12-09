from django.db import models
from edge.models import Fragment


class Operation(models.Model):
    class Meta:
        app_label = "edge"

    # numerical values are stored in database, do not change unless you have a
    # plan for migrating data.
    RECOMBINATION = (1, 'Homologous Recombination')
    PCR_SEQ_VERIFICATION = (2, 'PCR Product Sequence Verification')
    CRISPR_CUT = (3, 'CRISPR Cut')

    type = models.IntegerField(choices=(RECOMBINATION, PCR_SEQ_VERIFICATION, CRISPR_CUT))
    notes = models.TextField(null=True, blank=True)
    params = models.TextField(null=True, blank=True)
    fragment = models.ForeignKey(Fragment)
