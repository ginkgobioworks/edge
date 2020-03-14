import json
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

    @property
    def feature_name(self):
        if self.feature_base_first != 1 or self.feature_base_last != self.feature.length:
            return '%s[%s:%s]' % (self.feature.name,
                                  self.feature_base_first,
                                  self.feature_base_last)
        else:
            return self.feature.name

    def __str__(self):
        s = []
        s.append(self.feature_name)
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
               annotations[-1].feature_base_last == cf.feature_base_first - 1 and\
               annotations[-1].base_last == fcl.base_first - 1:
                # merge annotation
                annotations[-1].base_last = fcl.base_last
                annotations[-1].feature_base_last = cf.feature_base_last
            else:
                annotations.append(Annotation(base_first=fcl.base_first,
                                              base_last=fcl.base_last,
                                              chunk_feature=cf,
                                              fragment=fcl.fragment))
        return sorted(annotations, key=lambda a: a.base_first)


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
            try:
                self.id = klass.objects.select_for_update().order_by(
                    '-id').values('id')[0]['id'] + 1
            except IndexError:
                self.id = 1
        return super(BigIntPrimaryModel, self).save(*args, **kwargs)


class Chunk(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    initial_fragment = models.ForeignKey('Fragment', on_delete=models.PROTECT)
    sequence = models.TextField(null=True)

    def reload(self):
        return Chunk.objects.get(pk=self.pk)


class Edge(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    from_chunk = models.ForeignKey(Chunk, related_name='out_edges', on_delete=models.PROTECT)
    fragment = models.ForeignKey('Fragment', on_delete=models.PROTECT)
    # can be null, so we can supersede an edge from a child fragment
    to_chunk = models.ForeignKey(Chunk, null=True, related_name='in_edges',
                                 on_delete=models.PROTECT)


class Feature(models.Model):
    class Meta:
        app_label = "edge"

    name = models.CharField(max_length=100)
    type = models.CharField(max_length=100)
    strand = models.IntegerField(null=True)
    length = models.IntegerField()
    operation = models.ForeignKey('Operation', null=True, on_delete=models.CASCADE)
    _qualifiers = models.TextField(null=True, db_column='qualifiers')

    def set_qualifiers(self, qualifiers):
        self._qualifiers = json.dumps(qualifiers)

    @property
    def qualifiers(self):
        if self._qualifiers is not None:
            return json.loads(self._qualifiers)
        else:
            return None


class Chunk_Feature_Manager(models.Manager):
    def get_query_set(self):
        return super(Chunk_Feature_Manager, self).get_query_set().select_related('chunk', 'feature')


class Chunk_Feature(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"

    objects = Chunk_Feature_Manager()
    chunk = models.ForeignKey(Chunk, on_delete=models.PROTECT)
    feature = models.ForeignKey(Feature, on_delete=models.PROTECT)
    feature_base_first = models.IntegerField()
    feature_base_last = models.IntegerField()


class Fragment_Chunk_Location(BigIntPrimaryModel):
    class Meta:
        app_label = "edge"
        unique_together = (('fragment', 'chunk'),)
        index_together = (('fragment', 'base_last'), ('fragment', 'base_first'))

    fragment = models.ForeignKey('Fragment', on_delete=models.PROTECT)
    chunk = models.ForeignKey(Chunk, on_delete=models.PROTECT)
    base_first = models.IntegerField()
    base_last = models.IntegerField()

    def __str__(self):
        return '%s: %s-%s' % (self.fragment.name, self.location[0], self.location[1])

    @property
    def next_chunk(self):
        for fcl in self.fragment.fragment_chunk_location_set\
                                .select_related('chunk')\
                                .filter(base_first=self.base_last + 1):
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
        fc = Fragment_Chunk_Location.objects.get(fragment=self.fragment, base_last=loc[0] - 1)
        return fc

    def annotations(self):
        return [Annotation(base_first=self.base_first, base_last=self.base_last, chunk_feature=cf)
                for cf in self.chunk.chunk_feature_set.all()]

    @staticmethod
    @transaction.atomic()  # this method has to be in a transaction, see below
    def bulk_create(entries):
        # IMPORTANT: this entire method must be within a transaction, so we can
        # compute and assign unique BigIntPrimary IDs

        cur_id = 1
        try:
            cur_id = Fragment_Chunk_Location.objects\
                                            .select_for_update()\
                                            .order_by('-id')\
                                            .values('id')[0]['id'] + 1
        except IndexError:
            pass

        for entry in entries:
            entry.id = cur_id
            cur_id += 1
        Fragment_Chunk_Location.objects.bulk_create(entries)
