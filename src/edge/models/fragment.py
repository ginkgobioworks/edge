from django.db import models
from django.db import transaction


class Annotation(object):
    """
    Can contain multiple Chunk_Feature objects merged together.
    """

    def __init__(self, base_first, base_last, chunk_feature, fragment=None):
        self.base_first = base_first
        self.base_last = base_last
        self.feature = chunk_feature.feature
        self.feature_base_first = chunk_feature.feature_base_first
        self.feature_base_last = chunk_feature.feature_base_last
        self.fragment = fragment

    def __unicode__(self):
        s = []
        if self.feature_base_first != 1 or self.feature_base_last != self.feature.length:
            s.append('%s (%s-%s)' % (self.feature.name,
                                     self.feature_base_first,
                                     self.feature_base_last))
        else:
            s.append(self.feature.name)
        s.append(self.feature.type)
        if self.feature.strand == 1:
            s.append('+')
        else:
            s.append('-')
        return ', '.join(s)

    @staticmethod
    def from_chunk_feature_and_location_array(chunk_feature_locs):
        """
        chunk_feature_locs: an array of (Chunk_Feature, Fragment_Chunk_Location) tuples.
        """

        chunk_feature_locs = sorted(chunk_feature_locs,
                                    key=lambda t: (t[0].feature.id, t[1].base_first))

        annotations = []
        for cf, fcl in chunk_feature_locs:
            if len(annotations) > 0 and\
               annotations[-1].feature.id == cf.feature_id and\
               annotations[-1].feature_base_last == cf.feature_base_first-1 and\
               annotations[-1].base_last == fcl.base_first-1:
                # merge annotation
                annotations[-1].base_last = fcl.base_last
                annotations[-1].feature_base_last = cf.feature_base_last
            else:
                annotations.append(Annotation(base_first=fcl.base_first,
                                              base_last=fcl.base_last,
                                              chunk_feature=cf,
                                              fragment=fcl.fragment))
        return annotations


class BigIntPrimaryModel(models.Model):
    class Meta:
        app_label = "edge"
        abstract = True

    id = models.BigIntegerField(primary_key=True)

    @transaction.atomic()
    def save(self, *args, **kwargs):
        # mimic auto_increment
        if self.id is None:
            klass = type(self)
            if klass.objects.count() > 0:
                self.id = klass.objects.select_for_update().order_by('-id').values('id')[0]['id']+1
            else:
                self.id = 1
        return super(BigIntPrimaryModel, self).save(*args, **kwargs)


class Fragment(models.Model):
    class Meta:
        app_label = "edge"

    circular = models.BooleanField()
    name = models.CharField(max_length=256)
    parent = models.ForeignKey('self', null=True)
    start_chunk = models.ForeignKey('Chunk', null=True)

    @property
    def has_location_index(self):
        return self.fragment_chunk_location_set.count() > 0

    @property
    def length(self):
        q = self.fragment_chunk_location_set.order_by('-base_last')[:1]
        q = list(q)
        if len(q) == 0:
            return 0
        else:
            return q[0].base_last

    def predecessors(self):
        pred = [self]
        f = self.parent
        while f is not None:
            pred.append(f)
            f = f.parent
        return pred

    def predecessor_priorities(self):
        return {f.id: i for i, f in enumerate(self.predecessors())}

    def next_chunk(self, chunk):
        # check inheritance hierarchy to figure out which edge to use.
        if chunk.out_edges.count() == 0:
            return none
        else:
            # sort edges by predecessor level if more than one edge
            if chunk.out_edges.count() > 1:
                pp = self.predecessor_priorities()

                def sorter_f(e):
                    return pp[e.fragment_id] if e.fragment_id in pp else len(pp)

                out_edges = sorted(list(chunk.out_edges.all()), key=sorter_f)

            else:
                out_edges = chunk.out_edges.all()

            return out_edges[0].to_chunk

    def fragment_chunk(self, chunk):
        return self.fragment_chunk_location_set.filter(chunk=chunk)[0]

    def chunks(self, force_walk=False):
        if force_walk is False and self.has_location_index:
            q = self.fragment_chunk_location_set.select_related('chunk').order_by('base_first')
            for fcl in q:
                yield fcl.chunk
        else:
            chunk = self.start_chunk
            while chunk is not None:
                yield chunk
                chunk = self.next_chunk(chunk)

    @transaction.atomic()
    def index_fragment_chunk_locations(self):
        # remove old index
        self.fragment_chunk_location_set.all().delete()

        # go through each chunk and add new index
        i = 1
        for chunk in self.chunks(force_walk=True):
            if len(chunk.sequence) > 0:
                self.fragment_chunk_location_set.create(
                    chunk=chunk, base_first=i, base_last=i+len(chunk.sequence)-1
                )
                i += len(chunk.sequence)

    def get_sequence(self, bp_lo=None, bp_hi=None):
        q = self.fragment_chunk_location_set.select_related('chunk')
        if bp_lo is not None:
            q = q.filter(base_last__gte=bp_lo)
        if bp_hi is not None:
            q = q.filter(base_first__lte=bp_hi)
        q = q.order_by('base_first')

        sequence = []
        last_chunk_base_last = None

        for fcl in q:
            s = fcl.chunk.sequence
            if last_chunk_base_last is not None and fcl.base_first != last_chunk_base_last+1:
                raise Exception('Fragment chunk location table missing chunks before %s'
                                % (fcl.base_first,))
            if bp_lo is not None and fcl.base_first < bp_lo:
                s = s[bp_lo-fcl.base_first:]
            if bp_hi is not None and fcl.base_last > bp_hi:
                s = s[:bp_hi-fcl.base_last]
            sequence.append(s)
            last_chunk_base_last = fcl.base_last

        return ''.join(sequence)

    @property
    def sequence(self):
        return self.get_sequence()

    def annotations(self, bp_lo=None, bp_hi=None):
        q = self.fragment_chunk_location_set.select_related('chunk')
        if bp_lo is not None:
            q = q.filter(base_last__gte=bp_lo)
        if bp_hi is not None:
            q = q.filter(base_first__lte=bp_hi)
        q = q.order_by('base_first')

        chunk_features = []
        for fcl in q:
            feature_bps = [(f, fcl) for f in list(fcl.chunk.chunk_feature_set.all())]
            chunk_features.extend(feature_bps)

        return Annotation.from_chunk_feature_and_location_array(chunk_features)

    def update(self, name):
        from edge.fragment_writer import Fragment_Updater
        new_fragment = Fragment_Updater(
            name=name, circular=self.circular, parent=self, start_chunk=self.start_chunk
        )
        new_fragment.save()

        if self.has_location_index:
            # copy over location index
            for fc in self.fragment_chunk_location_set.all():
                new_fragment.fragment_chunk_location_set.create(
                    chunk=fc.chunk,
                    base_first=fc.base_first,
                    base_last=fc.base_last
                )
        else:
            new_fragment.index_fragment_chunk_locations()

        return new_fragment

    def annotate(self):
        if not self.has_location_index:
            self.index_fragment_chunk_locations()
        from edge.fragment_writer import Fragment_Annotator
        return Fragment_Annotator.objects.get(pk=self.pk)

    @staticmethod
    def non_genomic_fragments():
        fragments = list(Fragment.objects.filter(genome_fragment__id__isnull=True))
        return fragments

    # After some trials, we found that on MySQL, initial chunk size set to 20K
    # produced the best import time for the E. coli genome. Setting an initial
    # chunk size significantly improves import time: when dividing up smaller
    # chunks, you copy/slice less sequence data.
    @staticmethod
    def create_with_sequence(name, sequence, circular=False, initial_chunk_size=20000):
        from edge.fragment_writer import Fragment_Updater
        new_fragment = Fragment_Updater(name=name, circular=circular, parent=None, start_chunk=None)
        new_fragment.save()
        if initial_chunk_size is None or initial_chunk_size == 0:
            new_fragment.insert_bases(None, sequence)
        else:
            for i in range(0, len(sequence), initial_chunk_size):
                new_fragment.insert_bases(None, sequence[i:i+initial_chunk_size])
        return Fragment.objects.get(pk=new_fragment.pk)


class Chunk(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    initial_fragment = models.ForeignKey(Fragment)
    sequence = models.TextField(null=True)

    def reload(self):
        return Chunk.objects.get(pk=self.pk)


class Edge(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    from_chunk = models.ForeignKey(Chunk, related_name='out_edges')
    fragment = models.ForeignKey(Fragment)
    # can be null, so we can supersede an edge from a child fragment
    to_chunk = models.ForeignKey(Chunk, null=True, related_name='in_edges')


class Feature(models.Model):
    class Meta:
        app_label = "edge"

    name = models.CharField(max_length=100)
    type = models.CharField(max_length=100)
    strand = models.IntegerField(null=True)
    length = models.IntegerField()


class Chunk_Feature_Manager(models.Manager):
    def get_query_set(self):
        return super(Chunk_Feature_Manager, self).get_query_set().select_related('chunk', 'feature')


class Chunk_Feature(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    objects = Chunk_Feature_Manager()
    chunk = models.ForeignKey(Chunk)
    feature = models.ForeignKey(Feature)
    feature_base_first = models.IntegerField()
    feature_base_last = models.IntegerField()


class Fragment_Chunk_Location(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"
        unique_together = (('fragment', 'chunk'),)
        index_together = (('fragment', 'base_last'), ('fragment', 'base_first'))

    fragment = models.ForeignKey(Fragment)
    chunk = models.ForeignKey(Chunk)
    base_first = models.IntegerField()
    base_last = models.IntegerField()

    def __unicode__(self):
        return '%s: %s-%s' % (self.fragment.name, self.location[0], self.location[1])

    @property
    def next_chunk(self):
        for fcl in self.fragment.fragment_chunk_location_set\
                                .select_related('chunk')\
                                .filter(base_first=self.base_last+1):
            return fcl.chunk
        return None

    @property
    def location(self):
        return (self.base_first, self.base_last)

    @property
    def prev_fragment_chunk(self):
        loc = self.location
        if loc[0] == 1:
            return None
        fc = Fragment_Chunk_Location.objects.get(fragment=self.fragment, base_last=loc[0]-1)
        return fc

    def annotations(self):
        return [Annotation(base_first=self.base_first, base_last=self.base_last, chunk_feature=cf)
                for cf in self.chunk.chunk_feature_set.all()]
