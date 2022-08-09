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
    parent = models.ForeignKey(
        "self", null=True, on_delete=models.PROTECT, related_name="children"
    )
    notes = models.TextField(null=True, blank=True)
    fragments = models.ManyToManyField(Fragment, through="Genome_Fragment")
    created_on = models.DateTimeField("Created", auto_now_add=True, null=True)
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
    def import_gff(name, gff_fasta_fn, dirn='.'):
        genome = Genome.create(name)
        GFFImporter(genome, gff_fasta_fn).do_import(dirn=dirn)
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
        for f in sorted(self.fragments.all(), key=lambda f: f.id):
            f.indexed_fragment()
        # casting to Indexed_Genome
        return Indexed_Genome.objects.get(pk=self.pk)

    def root(self):
        if self.parent is None:
            return self
        return self.parent.root()

    def lock(self):
        """
        Do a select for update, which, when executed inside a transaction, places
        a lock on the genome.
        """

        root = self.root()

        genomes = Genome.objects.select_for_update().filter(pk=root.id)
        # Lock only happens when queryset is evaluated, therefore need to do at least genomes[0]
        genome = genomes[0]
        print(f"Lock genome {genome.id}")
        return self


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
            for fcl in Fragment_Chunk_Location.objects.filter(
                chunk=cf.chunk, fragment__genome=self
            ):
                cf_fcl.append((cf, fcl))
        annotations = Annotation.from_chunk_feature_and_location_array(cf_fcl)

        by_f = {}
        for annotation in annotations:
            if annotation.fragment.id not in by_f:
                by_f[annotation.fragment.id] = []
            by_f[annotation.fragment.id].append(annotation)
        return by_f

    def find_annotation_by_feature(self, feature):
        q = Chunk_Feature.objects.filter(
            chunk__fragment_chunk_location__fragment__genome=self, feature=feature
        )
        return self.__annotations_from_chunk_features(list(q))

    def find_annotation_by_name(self, name):
        q = Chunk_Feature.objects.filter(
            chunk__fragment_chunk_location__fragment__genome=self, feature__name=name
        )
        return self.__annotations_from_chunk_features(list(q))

    def find_annotation_by_qualifier(self, name, fields=None):
        fields = None if fields is None else [f.lower() for f in fields]

        q = Chunk_Feature.objects.filter(
            chunk__fragment_chunk_location__fragment__genome=self,
            feature___qualifiers__icontains=name,
        )
        chunk_features = []
        for cf in q:
            qualifiers = cf.feature.qualifiers
            for k, v in qualifiers.items():
                if fields is None or k.lower() in fields:
                    if type(v) in (bytes, str):
                        v = v.split(",")
                    v = [s.lower() for s in v]
                    if name.lower() in v:
                        chunk_features.append(cf)

        return self.__annotations_from_chunk_features(chunk_features)

    def changes(self):
        if self.parent is None:
            return []

        edges = Edge.objects.filter(
            fragment__genome_fragment__genome=self,
            fragment__genome_fragment__inherited=False,
        )

        fcs = []
        added = []
        for edge in edges:
            if edge.from_chunk:
                if edge.from_chunk.id not in added:
                    fcs.append(
                        edge.fragment.indexed_fragment().fragment_chunk(edge.from_chunk)
                    )
                    added.append(edge.from_chunk.id)
            if edge.to_chunk:
                if edge.to_chunk.id not in added:
                    fcs.append(
                        edge.fragment.indexed_fragment().fragment_chunk(edge.to_chunk)
                    )
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

    def get_coordinate_diff_from_parent_genome(self, p_genome_id):
        # Get fragments contained by child and not by parent
        p_fragment_ids = Fragment.objects.filter(genome=p_genome_id).values_list('id', flat=True)
        child_only_fragments = self.fragments.filter().exclude(genome=p_genome_id)

        # Relate child fragments to parent fragments
        child_to_parent_fs = {}
        for f in child_only_fragments:
            curr_f = f
            while curr_f is not None and curr_f.id not in p_fragment_ids:
                curr_f = curr_f.parent
            if curr_f is not None:
                child_to_parent_fs[f.indexed_fragment()] = curr_f.indexed_fragment()
        if len(child_to_parent_fs) == 0:
            raise Exception("Genome input is not a parent by lineage")

        # Get start to end chunks by fragment of diff regions between parent and child
        fragment_to_starts_and_ends = {}
        for child_f, parent_f in child_to_parent_fs.items():
            parent_chunks, child_chunks = list(parent_f.chunks()), list(child_f.chunks())
            starts_and_ends, current_start = {}, None
            for child_i, c in enumerate(child_chunks):
                if c in parent_chunks:
                    parent_i = parent_chunks.index(c)
                    prev_parent_chunk = parent_chunks[parent_i - 1] \
                        if parent_i > 0 else None
                    next_parent_chunk = parent_chunks[parent_i + 1] \
                        if parent_i < (len(parent_chunks) - 1) else None
                    prev_child_chunk = child_chunks[child_i - 1] \
                        if child_i > 0 else None
                    next_child_chunk = child_chunks[child_i + 1] \
                        if child_i < (len(child_chunks) - 1) else None

                    if prev_parent_chunk != prev_child_chunk:
                        if current_start is not None:
                            starts_and_ends[current_start] = c.id
                            current_start = None
                        if len(starts_and_ends) == 0:
                            starts_and_ends[None] = c.id

                    if next_parent_chunk != next_child_chunk:
                        current_start = c.id

            if current_start is not None:
                starts_and_ends[current_start] = None
            fragment_to_starts_and_ends[child_f] = starts_and_ends

        # Store dictionaries of changes with ID, start, and end for parent and child fragments
        regions = []
        for child_f, starts_and_ends in fragment_to_starts_and_ends.items():
            parent_f = child_to_parent_fs[child_f]
            for start_id, end_id in starts_and_ends.items():
                regions.append({
                    "parent_fragment_name": parent_f.name,

                    "parent_fragment_id": parent_f.id,

                    "parent_starts_at": 1 if start_id is None
                    else parent_f.fragment_chunk(start_id).base_last + 1,

                    "parent_ends_before": parent_f.length + 1 if end_id is None
                    else parent_f.fragment_chunk(end_id).base_first,

                    "child_fragment_name": child_f.name,

                    "child_fragment_id": child_f.id,

                    "child_starts_at": 1 if start_id is None
                    else child_f.fragment_chunk(start_id).base_last + 1,

                    "child_ends_before": child_f.length + 1 if end_id is None
                    else child_f.fragment_chunk(end_id).base_first,
                })

        return regions


class Genome_Fragment(models.Model):
    class Meta:
        app_label = "edge"

    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    fragment = models.ForeignKey(Fragment, on_delete=models.CASCADE)
    inherited = models.BooleanField()
