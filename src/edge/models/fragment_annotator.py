from django.db import transaction
from django.db.models import F
from edge.models.chunk import *


class Fragment_Writer:
    """
    Mixin that includes helpers for updating a fragment.
    """

    def _add_feature(self, name, type, length, strand):
        if strand not in (1, -1, None):
            raise Exception('Strand must be 1, -1, or None')
        f = Feature(name=name, type=type, length=length, strand=strand)
        f.save()
        return f

    def _annotate_chunk(self, chunk, feature, feature_base_first, feature_base_last):
        Chunk_Feature(chunk=chunk, feature=feature,
                      feature_base_first=feature_base_first,
                      feature_base_last=feature_base_last).save()

    def _add_chunk(self, sequence, fragment):
        c = Chunk(sequence=sequence, initial_fragment=fragment)
        c.save()
        return c

    def _reset_chunk_sequence(self, chunk, sequence):
        chunk.sequence = sequence
        chunk.save()
        Chunk_Feature.objects.filter(chunk=chunk).delete()
        # if chunk is the start chunk, update the object reference
        if self.start_chunk.id == chunk.id:
            self.start_chunk = chunk
            self.save()

    def _assert_not_linked_to(self, chunk):
        if Edge.objects.filter(to_chunk=chunk, fragment=self).count() > 0:
            raise Exception('Fragment %s already linked to chunk %s' % (self.id, chunk.id))

    def _add_edges(self, chunk, *unsaved_edges):
        existing_edges = list(Edge.objects.filter(from_chunk=chunk))

        # get list of fragment IDs in new edges
        new_fragment_ids = [e.fragment_id for e in unsaved_edges]

        # remove old edges with same fragment ids as new edges
        for e in existing_edges:
            if e.fragment_id in new_fragment_ids:
                e.delete()

        # add new edges
        for edge in unsaved_edges:
            edge.save()

    def _split_annotations(self, annotations, bps_to_split, split1, split2):
        for a in annotations:
            a1 = (a.feature, a.feature_base_first, a.feature_base_first+bps_to_split-1)
            a2 = (a.feature, a.feature_base_first+bps_to_split, a.feature_base_last)
            self._annotate_chunk(split1, *a1)
            self._annotate_chunk(split2, *a2)

    def _split_chunk(self, chunk, s1, s2):
        from edge.models.fragment import Fragment_Chunk_Location

        # splitted chunk should be "created" by the fragment that created the
        # original chunk
        split2 = self._add_chunk(s2, chunk.initial_fragment)
        self._reset_chunk_sequence(chunk, s1)

        # move all edges from original chunk to second chunk
        Edge.objects.filter(from_chunk=chunk).update(from_chunk=split2)
        # verify chunk edges been updated
        assert chunk.out_edges.count() == 0

        # add edge from chunk to split2
        self._add_edges(chunk, Edge(from_chunk=chunk,
                                    fragment=chunk.initial_fragment, to_chunk=split2))

        # add location for new chunk - note that the following code has to
        # update the location index for each of the fragment that uses the
        # chunk. there may be a lot of fragments, potentially a big scalability
        # problem.
        entries = []
        for fcl in chunk.fragment_chunk_location_set.all():
            entries.append(Fragment_Chunk_Location(fragment_id=fcl.fragment_id,
                                                   chunk_id=split2.id,
                                                   base_first=fcl.base_first+len(s1),
                                                   base_last=fcl.base_last))
        Fragment_Chunk_Location.bulk_create(entries)

        # adjust chunk location index for existing chunk
        chunk.fragment_chunk_location_set.update(base_last=F('base_first')+len(s1)-1)
        return split2

    def _find_chunk_prev_next(self, before_base1):
        prev_chunk = None
        next_chunk = None
        chunk = None
        bases_visited = None

        if before_base1 is not None:
            q = self.fragment_chunk_location_set.filter(base_first__lte=before_base1,
                                                        base_last__gte=before_base1)
            q = list(q)
            if len(q) > 0:
                fc = q[0]
                chunk = fc.chunk
                bases_visited = fc.base_last
                next_chunk = fc.next_chunk

                # find prev chunk
                if fc.base_first == 1:
                    prev_chunk = None
                else:
                    prev_chunk = fc.prev_fragment_chunk.chunk

        if chunk is None:  # after all sequence ended, need last chunk and total bases
            total_bases = self.length
            if total_bases > 0:
                prev_chunk = self.fragment_chunk_location_set.get(base_last=total_bases).chunk
            bases_visited = total_bases

        return prev_chunk, chunk, next_chunk, bases_visited

    def _find_and_split_before(self, before_base1):

        if before_base1 is not None and before_base1 <= 0:
            raise Exception('chunk index should be 1-based')

        chunk = None
        chunk_id = None
        prev_chunk = None
        next_chunk = None
        bases_visited = 0

        prev_chunk, chunk, next_chunk, bases_visited = self._find_chunk_prev_next(before_base1)

        if chunk:
            chunk_id = chunk.id

        # found the bp we are looking for
        if before_base1 is not None and bases_visited >= before_base1:
            chunk_len = len(chunk.sequence)

            # can avoid splitting if first bp in this chunk is before_base1
            if bases_visited-chunk_len+1 == before_base1:
                #print 'no need to split before %s' % (before_base1,)
                return prev_chunk, chunk

            # otherwise, have to split the chunk
            first_bp_in_chunk = bases_visited-chunk_len+1
            bps_to_split = before_base1-first_bp_in_chunk
            s1 = chunk.sequence[0:bps_to_split]
            s2 = chunk.sequence[bps_to_split:]

            # save original annotations, which will be trashed in _split_chunk (which
            # calls _reset_chunk_sequence)
            original_annotations = self.fragment_chunk(chunk).annotations()
            # split chunk
            split2 = self._split_chunk(chunk, s1, s2)
            chunk = chunk.reload()

            # split up annotations as well
            self._split_annotations(original_annotations, bps_to_split, chunk, split2)
            return chunk, split2

        else:  # found end of the fragment
            return prev_chunk, None


class Fragment_Annotator:
    """
    Mixin for annotating a fragment.
    """

    @transaction.atomic()
    def annotate(self, first_base1, last_base1, name, type, strand):
        if self.circular and last_base1 < first_base1:
            # has to figure out the total length from last chunk
            length = self.length-first_base1+1+last_base1
        else:
            length = last_base1-first_base1+1
            if length <= 0:
                raise Exception('Annotation must have length one or more')

        prev_chunk, annotation_start = self._find_and_split_before(first_base1)
        annotation_end, next_chunk = self._find_and_split_before(last_base1+1)

        # did two splits, so must reload annotation_start in case that got splitted
        annotation_start = annotation_start.reload()

        new_feature = self._add_feature(name, type, length, strand)

        # now, starting with chunk annotation_start, walk through chunks until
        # we hit annotation_end, and add annotation for each chunk
        chunk = annotation_start
        a_i = 1
        while True:
            fc = self.fragment_chunk(chunk)
            self._annotate_chunk(chunk, new_feature, a_i, a_i+len(chunk.sequence)-1)
            a_i += len(chunk.sequence)
            if chunk.id == annotation_end.id:
                break
            chunk = fc.next_chunk
            if chunk is None:
                chunk = self.start_chunk
