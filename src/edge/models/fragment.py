from django.db import models


class Annotation(object):
    """
    Can contain multiple Chunk_Feature objects merged together.
    """

    def __init__(self, base_first, base_last, chunk_feature):
        self.base_first = base_first
        self.base_last = base_last
        self.feature = chunk_feature.feature
        self.feature_base_first = chunk_feature.feature_base_first
        self.feature_base_last = chunk_feature.feature_base_last

    def __unicode__(self):
        s = []
        if self.feature_base_first != 1 or self.feature_base_last != self.feature.length:
            s.append('%s (%s-%s)' % (self.feature.name, self.feature_base_first, self.feature_base_last))
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

        chunk_feature_locs = sorted(chunk_feature_locs, key=lambda t: (t[0].feature.id, t[1].base_first))

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
                                              chunk_feature=cf))
        return annotations


class Fragment(models.Model):
    class Meta:
        app_label = "edge"

    circular = models.BooleanField()
    name = models.CharField(max_length=256)
    parent = models.ForeignKey('self', null=True)
    start_chunk = models.ForeignKey('Chunk', null=True)

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

    def fragment_chunk(self, chunk):
        return self.fragment_chunk_location_set.filter(chunk=chunk)[0]

    def chunks(self):
        q = self.fragment_chunk_location_set.select_related('chunk').order_by('base_first')
        for fcl in q:
          yield fcl.chunk

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

        # copy over location index
        for fc in self.fragment_chunk_location_set.all():
           new_fragment.fragment_chunk_location_set.create(
               chunk=fc.chunk,
               base_first=fc.base_first,
               base_last=fc.base_last
           )

        return new_fragment

    def annotate(self):
        from edge.fragment_writer import Fragment_Annotator
        return Fragment_Annotator.objects.get(pk=self.pk)

    @staticmethod
    def non_genomic_fragments():
        fragments = list(Fragment.objects.filter(genome_fragment__id__isnull=True))
        return fragments

    @staticmethod
    def create_with_sequence(name, sequence, circular=False):
        from edge.fragment_writer import Fragment_Updater
        new_fragment = Fragment_Updater(name=name, circular=circular, parent=None, start_chunk=None)
        new_fragment.save()
        new_fragment.insert_bases(None, sequence)
        new_fragment.save()
        return Fragment.objects.get(pk=new_fragment.pk)


class Chunk(models.Model):
    class Meta:
        app_label = "edge"

    id = models.BigIntegerField(primary_key=True)
    initial_fragment = models.ForeignKey(Fragment)
    sequence = models.TextField(null=True)

    def reload(self):
        return Chunk.objects.get(pk=self.pk)

    def save(self, *args, **kwargs):
        # mimic auto_increment
        if self.id is None:
            if Chunk.objects.count() > 0:
                self.id = Chunk.objects.order_by('-id')[0].id+1
            else:
                self.id = 1
        return super(Chunk, self).save(*args, **kwargs)


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


class Chunk_Feature(models.Model):
    class Meta:
        app_label = "edge"

    objects = Chunk_Feature_Manager()
    id = models.BigIntegerField(primary_key=True)
    chunk = models.ForeignKey(Chunk)
    feature = models.ForeignKey(Feature)
    feature_base_first = models.IntegerField()
    feature_base_last = models.IntegerField()

    def save(self, *args, **kwargs):
        # mimic auto_increment
        if self.id is None:
            if Chunk_Feature.objects.count() > 0:
                self.id = Chunk_Feature.objects.order_by('-id')[0].id+1
            else:
                self.id = 1
        return super(Chunk_Feature, self).save(*args, **kwargs)


class Fragment_Chunk_Location(models.Model):
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
    def sequence(self):
        return self.chunk.sequence

    @property
    def initial_fragment(self):
        return self.chunk.initial_fragment

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
