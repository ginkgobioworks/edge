from django.db import models

class _Genome(models.Model):
    name = models.CharField(max_length=256)
    parent = models.ForeignKey('self', null=True)
    notes = models.TextField(null=True)
    fragments = models.ManyToManyField(_Fragment, through='_Genome_Fragment')


class _Genome_Fragment(models.Model):
    genome = models.ForeignKey(_Genome)
    fragment = models.ForeignKey(_Fragment)
    inherited = models.BooleanField()
    

