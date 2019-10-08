from django.utils import timezone
from django.db import models
from django.db.models import Q

from edge.models.chunk import (
    Annotation,
    Chunk_Feature,
    Fragment_Chunk_Location,
)
from edge.models.fragment_writer import Fragment_Writer
from edge.models.fragment_annotator import Fragment_Annotator
from edge.models.fragment_updater import Fragment_Updater


class Fragment(models.Model):
    class Meta:
        app_label = "edge"

    circular = models.BooleanField()
    name = models.CharField(max_length=256)
    parent = models.ForeignKey('self', null=True, on_delete=models.PROTECT)
    start_chunk = models.ForeignKey('Chunk', null=True, on_delete=models.PROTECT)
    est_length = models.IntegerField('Estimated length', null=True, blank=True)
    created_on = models.DateTimeField('Created', auto_now_add=True, null=True)
    active = models.BooleanField(default=True)

    @staticmethod
    def user_defined_fragments(q=None, f=None, l=None):  # noqa: E741
        qs = Fragment.objects.filter(genome_fragment__id__isnull=True, active=True)
        if q is not None:
            qs = qs.filter(q)
        f = 0 if f is None else f
        if l is None:
            qs = qs[f:]
        else:
            qs = qs[f:l]
        fragments = list(qs)
        return fragments

    # After some trials, we found that on MySQL, initial chunk size set to 20K
    # produced the best import time for the E. coli genome. Setting an initial
    # chunk size significantly improves import time: when dividing up smaller
    # chunks, you copy/slice less sequence data.
    @staticmethod
    def create_with_sequence(name, sequence, circular=False, initial_chunk_size=20000):
        new_fragment = Fragment(name=name, circular=circular, parent=None, start_chunk=None)
        new_fragment.save()
        new_fragment = new_fragment.indexed_fragment()
        if initial_chunk_size is None or initial_chunk_size == 0:
            new_fragment.insert_bases(None, sequence)
        else:
            for i in range(0, len(sequence), initial_chunk_size):
                new_fragment.insert_bases(None, sequence[i:i + initial_chunk_size])
        return new_fragment

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
        """
        Finds successor of the specified chunk. Uses fragment inheritance
        hierarchy to pick an edge from the outward edges of the specified
        chunk.
        """

        if chunk.out_edges.count() == 0:
            return None
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

    def index_fragment_chunk_locations(self):
        # remove old index
        self.fragment_chunk_location_set.all().delete()

        # go through each chunk and add new index
        i = 1
        entries = []
        for chunk in self.chunks_by_walking():
            if len(chunk.sequence) > 0:
                entries.append(Fragment_Chunk_Location(fragment_id=self.id,
                                                       chunk_id=chunk.id,
                                                       base_first=i,
                                                       base_last=i + len(chunk.sequence) - 1))

                i += len(chunk.sequence)
        Fragment_Chunk_Location.bulk_create(entries)

        indexed = Indexed_Fragment.objects.get(pk=self.id)
        if indexed.length != self.est_length:
            self.est_length = indexed.length
            self.save()
            indexed.est_length = self.est_length

        try:
            index = self.fragment_index
        except BaseException:
            index = Fragment_Index(fragment=self)
        index.fresh = True
        index.updated_on = timezone.now()
        index.save()

        return indexed

    def chunks_by_walking(self):
        chunk = self.start_chunk
        while chunk is not None:
            yield chunk
            chunk = self.next_chunk(chunk)

    @property
    def has_location_index(self):
        if self.fragment_chunk_location_set.count() > 0:
            try:
                index = self.fragment_index
            except BaseException:
                return False
            return index.fresh
        return False

    def indexed_fragment(self):
        if not self.has_location_index:
            return self.index_fragment_chunk_locations()
        # casting to Indexed_Fragment
        return Indexed_Fragment.objects.get(pk=self.id)


class Fragment_Index(models.Model):
    class Meta:
        app_label = "edge"

    fragment = models.OneToOneField(Fragment, on_delete=models.CASCADE)
    fresh = models.BooleanField()
    updated_on = models.DateTimeField('Updated', null=True)


class Indexed_Fragment(Fragment_Annotator, Fragment_Updater, Fragment_Writer, Fragment):
    """
    An Indexed_Fragment is a Fragment with chunk location index. You need chunk
    location index to efficiently find annotations and bp positions.

    Invariant: Fragment#prepare method always return an Indexed_Fragment with
    location index.
    """

    class Meta:
        app_label = "edge"
        proxy = True

    def chunks(self):
        q = self.fragment_chunk_location_set.select_related('chunk').order_by('base_first')
        for fcl in q:
            yield fcl.chunk

    def circ_bp(self, bp):
        if self.circular is True:
            return ((bp - 1) % self.length) + 1
        return bp

    @property
    def length(self):
        q = self.fragment_chunk_location_set.order_by('-base_last')[:1]
        q = list(q)
        if len(q) == 0:
            return 0
        else:
            return q[0].base_last

    def fragment_chunk(self, chunk):
        return self.fragment_chunk_location_set.filter(chunk=chunk)[0]

    def __get_linear_sequence(self, bp_lo=None, bp_hi=None):
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
            if last_chunk_base_last is not None and fcl.base_first != last_chunk_base_last + 1:
                raise Exception('Fragment chunk location table missing chunks before %s'
                                % (fcl.base_first,))
            if bp_lo is not None and fcl.base_first < bp_lo:
                s = s[bp_lo - fcl.base_first:]
            if bp_hi is not None and fcl.base_last > bp_hi:
                s = s[:bp_hi - fcl.base_last]
            sequence.append(s)
            last_chunk_base_last = fcl.base_last

        return ''.join(sequence)

    def get_sequence(self, bp_lo=None, bp_hi=None):
        if bp_lo is not None:
            bp_lo = self.circ_bp(bp_lo)
            bp_lo = 1 if bp_lo < 1 else bp_lo

        if bp_hi is not None:
            bp_hi = self.circ_bp(bp_hi)
            bp_hi = self.length if bp_hi > self.length else bp_hi

        if bp_lo is None or bp_hi is None or bp_lo <= bp_hi:
            return self.__get_linear_sequence(bp_lo=bp_lo, bp_hi=bp_hi)

        else:
            assert self.circular is True
            s1 = self.__get_linear_sequence(bp_lo=bp_lo)
            s2 = self.__get_linear_sequence(bp_hi=bp_hi)
            return s1 + s2

    @property
    def sequence(self):
        return self.get_sequence()

    def annotations(self, bp_lo=None, bp_hi=None):
        # hand-crafted query to fetch both Chunk_Feature and
        # Fragment_Chunk_Location instances

        # because chunk__fragment_chunk_location is a M2M relationship, we have
        # to specify all the filtering rules in one filter method call. if the
        # rules are spread across multiple method calls, Django will construct
        # a separate join to the fragment_chunk_location table for each filter
        # call.
        rules = [Q(chunk__fragment_chunk_location__fragment=self)]
        if bp_lo is not None:
            rules.append(Q(chunk__fragment_chunk_location__base_last__gte=bp_lo))
        if bp_hi is not None:
            rules.append(Q(chunk__fragment_chunk_location__base_first__lte=bp_hi))

        fcl_tb = Fragment_Chunk_Location._meta.db_table
        bf = [f.column for f in Fragment_Chunk_Location._meta.fields if f.name == 'base_first'][0]
        bl = [f.column for f in Fragment_Chunk_Location._meta.fields if f.name == 'base_last'][0]
        q = Chunk_Feature.objects.filter(*rules)\
                                 .extra(select=dict(fcl_base_first='%s.%s' % (fcl_tb, bf),
                                                    fcl_base_last='%s.%s' % (fcl_tb, bl)))

        # using fcl_base_first and fcl_base_last fields to create a fake,
        # unsaved Fragment_Chunk_Location object
        chunk_features = [(cf, Fragment_Chunk_Location(fragment=self, chunk=cf.chunk,
                                                       base_first=int(cf.fcl_base_first),
                                                       base_last=int(cf.fcl_base_last)))
                          for cf in q]

        chunk_features = sorted(chunk_features, key=lambda t: t[1].base_first)
        return Annotation.from_chunk_feature_and_location_array(chunk_features)

    def update(self, name):
        new_fragment = \
            Fragment(name=name, circular=self.circular, parent=self, start_chunk=self.start_chunk)
        new_fragment.save()

        # copy over location index
        entries = []
        for fc in self.fragment_chunk_location_set.all():
            entries.append(Fragment_Chunk_Location(fragment_id=new_fragment.id,
                                                   chunk_id=fc.chunk_id,
                                                   base_first=fc.base_first,
                                                   base_last=fc.base_last))
        Fragment_Chunk_Location.bulk_create(entries)
        Fragment_Index(fragment=new_fragment, fresh=True,
                       updated_on=self.fragment_index.updated_on).save()
        return new_fragment.indexed_fragment()
