from django.db import models
from edge.models.fragment import *


class Genome(models.Model):
    class Meta:
        app_label = "edge"

    name = models.CharField(max_length=256)
    parent = models.ForeignKey('self', null=True)
    notes = models.TextField(null=True)
    fragments = models.ManyToManyField(Fragment, through='Genome_Fragment')

    @staticmethod
    def create(name, notes=None):
        new_genome = Genome(name=name, notes=notes, parent=None)
        new_genome.save()
        return Genome_Updater.objects.get(pk=new_genome.pk)

    def edit(self):
        from edge.genome_writer import Genome_Updater
        return Genome_Updater.objects.get(pk=self.pk)

    def annotate(self):
        return self.edit()

    def update(self, name=None, notes=None):
        new_genome = Genome(name=name, notes=notes, parent=self)
        new_genome.save()
        for f in self.fragments.all():
          Genome_Fragment(genome=new_genome, fragment=f, inherited=True).save()
        return Genome_Updater.objects.get(pk=new_genome.pk)

    def find_annotation(self, name):
        q = Chunk_Feature.objects.filter(chunk__fragment_chunk_location__fragment__genome=self,
                                         feature__name=name)
        chunk_features = []
        for cf in q:
          # get fragment chunk location
          for fcl in Fragment_Chunk_Location.objects.filter(chunk=cf.chunk, fragment__genome=self):
            chunk_features.append(cf, fcl)
        return Annotation.from_chunk_feature_and_location_array(chunk_features)

    def changes(self):
        if self.parent is None:
            return []

        edges = Edge.objects.filter(
          fragment__genome_fragment__genome=self,
          fragment__genome_fragment__inherited=False
        )

        # XXX isn't this incredibly slow?
        # XXX remove duplicate chunks?
        chunks = []
        for edge in edges:
          if edge.from_chunk:
            chunks.append(edge.fragment.fragment_chunk(edge.from_chunk))
          if edge.to_chunk:
            chunks.append(edge.fragment.fragment_chunk(edge.to_chunk))

        return chunks

    def changed_locations_by_fragment(self):
        changes = {}

        for c in self.changes():
            if c.fragment.id not in changes:
                changes[c.fragment.id] = []
            v = changes[c.fragment.id]
            if len(v) == 0 or v[-1][1]+1 != c.location[0]:
                v.append([c.location[0], c.location[1]])
            else:
                v[-1][1] = c.location[1]

        changes = {self.get_fragment_by_id(f): v for f, v in changes.iteritems()}
        return changes


class Genome_Fragment(models.Model):
    class Meta:
        app_label = "edge"

    genome = models.ForeignKey(Genome)
    fragment = models.ForeignKey(Fragment)
    inherited = models.BooleanField()
