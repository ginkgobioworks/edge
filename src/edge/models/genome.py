from django.db import models

from edge.importer import GFFImporter
from edge.models.chunk import (
    Annotation,
    Chunk_Feature,
    Edge,
    Fragment_Chunk_Location,
)
from edge.models.fragment import Fragment
from edge.models.genome_updater import Genome_Updater


class Genome(Genome_Updater, models.Model):
    class Meta:
        app_label = "edge"

    name = models.TextField()
    parent = models.ForeignKey('self', null=True, on_delete=models.PROTECT,
                               related_name='children')
    notes = models.TextField(null=True, blank=True)
    fragments = models.ManyToManyField(Fragment, through='Genome_Fragment')
    created_on = models.DateTimeField('Created', auto_now_add=True, null=True)
    active = models.BooleanField(default=True)
    blastdb = models.TextField(null=True, blank=True)

    def __str__(self):
        return self.name

    @staticmethod
    def create(name, notes=None):
        new_genome = Genome(name=name, notes=notes, parent=None)
        new_genome.save()
        return new_genome

    @staticmethod
    def import_gff(name, gff_fasta_fn):
        genome = Genome.create(name)
        GFFImporter(genome, gff_fasta_fn).do_import()
        return genome

    def update(self, name=None, notes=None):
        name = self.name if name is None else name
        new_genome = Genome(name=name, notes=notes, parent=self)
        new_genome.save()
        for f in self.fragments.all():
            Genome_Fragment(genome=new_genome, fragment=f, inherited=True).save()
        return new_genome

    @property
    def has_location_index(self):
        for f in self.fragments.all():
            if not f.has_location_index:
                return False
        return True

    def indexed_genome(self):
        for f in self.fragments.all():
            f.indexed_fragment()
        # casting to Indexed_Genome
        return Indexed_Genome.objects.get(pk=self.pk)


class Indexed_Genome(Genome):
    """
    An Indexed_Genome is a Genome where each fragment has been indexed.
    """

    class Meta:
        app_label = "edge"
        proxy = True

    def __annotations_from_chunk_features(self, chunk_features):
        cf_fcl = []
        for cf in chunk_features:
            # get fragment chunk location
            for fcl in Fragment_Chunk_Location.objects.filter(chunk=cf.chunk,
                                                              fragment__genome=self):
                cf_fcl.append((cf, fcl))
        annotations = Annotation.from_chunk_feature_and_location_array(cf_fcl)

        by_f = {}
        for annotation in annotations:
            if annotation.fragment.id not in by_f:
                by_f[annotation.fragment.id] = []
            by_f[annotation.fragment.id].append(annotation)
        return by_f

    def find_annotation_by_feature(self, feature):
        q = Chunk_Feature.objects.filter(chunk__fragment_chunk_location__fragment__genome=self,
                                         feature=feature)
        return self.__annotations_from_chunk_features(list(q))

    def find_annotation_by_name(self, name):
        q = Chunk_Feature.objects.filter(chunk__fragment_chunk_location__fragment__genome=self,
                                         feature__name=name)
        return self.__annotations_from_chunk_features(list(q))

    def find_annotation_by_qualifier(self, name, fields=None):
        fields = None if fields is None else [f.lower() for f in fields]

        q = Chunk_Feature.objects.filter(chunk__fragment_chunk_location__fragment__genome=self,
                                         feature___qualifiers__icontains=name)
        chunk_features = []
        for cf in q:
            qualifiers = cf.feature.qualifiers
            for k, v in qualifiers.items():
                if fields is None or k.lower() in fields:
                    if type(v) in (bytes, str):
                        v = v.split(',')
                    v = [s.lower() for s in v]
                    if name.lower() in v:
                        chunk_features.append(cf)

        return self.__annotations_from_chunk_features(chunk_features)

    def changes(self):
        if self.parent is None:
            return []

        edges = Edge.objects.filter(
            fragment__genome_fragment__genome=self,
            fragment__genome_fragment__inherited=False
        )

        fcs = []
        added = []
        for edge in edges:
            if edge.from_chunk:
                if edge.from_chunk.id not in added:
                    fcs.append(edge.fragment.indexed_fragment().fragment_chunk(edge.from_chunk))
                    added.append(edge.from_chunk.id)
            if edge.to_chunk:
                if edge.to_chunk.id not in added:
                    fcs.append(edge.fragment.indexed_fragment().fragment_chunk(edge.to_chunk))
                    added.append(edge.to_chunk.id)
        return fcs

    def changed_locations_by_fragment(self):
        changes = {}

        for c in self.changes():
            if c.fragment.id not in changes:
                changes[c.fragment.id] = []
            v = changes[c.fragment.id]
            if len(v) == 0 or v[-1][1] + 1 != c.location[0]:
                v.append([c.location[0], c.location[1]])
            else:
                v[-1][1] = c.location[1]

        changes = {Fragment.objects.get(pk=f): v for f, v in changes.items()}
        return changes


class Genome_Fragment(models.Model):
    class Meta:
        app_label = "edge"

    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    fragment = models.ForeignKey(Fragment, on_delete=models.CASCADE)
    inherited = models.BooleanField()
