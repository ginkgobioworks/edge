from django.db import transaction
from django.db.models import F
from edge.models import *


class Fragment_Writer(Fragment):
    class Meta:
        app_label = "edge"
        proxy = True

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

    def _find_chunk_prev_next_by_location_index(self, before_base1):
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

        prev_chunk, chunk, next_chunk, bases_visited = \
            self._find_chunk_prev_next_by_location_index(before_base1)

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


class Fragment_Annotator(Fragment_Writer):
    class Meta:
        app_label = "edge"
        proxy = True

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


class Fragment_Updater(Fragment_Writer):
    class Meta:
        app_label = "edge"
        proxy = True

    def _append_to_fragment(self, prev_chunk, cur_fragment_length, sequence):
        # only use this if you are appending chunk to fragment while importing
        # a fragment
        new_chunk = self._add_chunk(sequence, self)
        if prev_chunk is not None:  # add chunks after prev_chunk_id
            self._add_edges(prev_chunk,
                            Edge(from_chunk=prev_chunk, fragment=self, to_chunk=new_chunk))
        else:
            self.start_chunk = new_chunk
            self.save()
        self._add_edges(new_chunk, Edge(from_chunk=new_chunk, fragment=self, to_chunk=None))
        self.fragment_chunk_location_set.create(
            chunk=new_chunk,
            base_first=cur_fragment_length+1,
            base_last=cur_fragment_length+1+len(sequence)-1
        )
        return new_chunk

    @transaction.atomic()
    def insert_bases(self, before_base1, sequence):
        # find chunks before and containing the insertion point
        prev_chunk, chunk = self._find_and_split_before(before_base1)

        # create new chunk
        new_chunk = self._add_chunk(sequence, self)

        if prev_chunk is not None:  # add chunks after prev_chunk_id
            self._add_edges(prev_chunk, Edge(from_chunk=prev_chunk,
                                             fragment=self, to_chunk=new_chunk))

        else:  # add chunks at start of fragment
            self.start_chunk = new_chunk
            self.save()

        # chunk may be None, but that's okay, we want to make sure this chunk
        # is the END and not going to be superseded by child fragment appending
        # more chunks!
        self._add_edges(new_chunk, Edge(from_chunk=new_chunk, fragment=self, to_chunk=chunk))

        # shift base_first and base_last for existing chunks
        if before_base1 is not None:
            self.fragment_chunk_location_set.filter(base_first__gte=before_base1)\
                                            .update(base_first=F('base_first')+len(sequence),
                                                    base_last=F('base_last')+len(sequence))

        # insert location for new chunk
        if before_base1 is not None:
            self.fragment_chunk_location_set.create(
                chunk=new_chunk,
                base_first=before_base1,
                base_last=before_base1+len(sequence)-1
            )
        else:
            fragment_length = self.length
            self.fragment_chunk_location_set.create(
                chunk=new_chunk,
                base_first=fragment_length+1,
                base_last=fragment_length+1+len(sequence)-1
            )

    @transaction.atomic()
    def remove_bases(self, before_base1, length):
        if length <= 0:
            raise Exception('Cannot remove less than one base pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        prev_chunk, dont_use_stale_after_second_split = self._find_and_split_before(before_base1)
        removal_end, next_chunk = self._find_and_split_before(before_base1+length)

        if prev_chunk is not None:  # remove chunks after prev_chunk_id
            if next_chunk:  # delete prior to end of fragment
                self._add_edges(prev_chunk,
                                Edge(from_chunk=prev_chunk, fragment=self, to_chunk=next_chunk))
            else:
                # delete all remaining chunks in fragment by adding an edge with None
                # as target, superseding any child that may have added chunks after
                # prev_chunk_id
                self._add_edges(prev_chunk,
                                Edge(from_chunk=prev_chunk, fragment=self, to_chunk=None))

        else:  # remove chunks at start of fragment
            if next_chunk is None:
                raise Exception('Cannot remove entire fragment')
            self.start_chunk = next_chunk
            self.save()

        # remove location for deleted chunks
        self.fragment_chunk_location_set.filter(base_first__gte=before_base1,
                                                base_first__lt=before_base1+length).delete()

        # shift base_first and base_last for existing chunks
        self.fragment_chunk_location_set.filter(base_first__gt=before_base1)\
                                        .update(base_first=F('base_first')-length,
                                                base_last=F('base_last')-length)

    @transaction.atomic()
    def replace_bases(self, before_base1, length_to_remove, sequence):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_bases(before_base1, sequence)

    @transaction.atomic()
    def insert_fragment(self, before_base1, fragment):

        # find chunks before and containing the insertion point
        prev_chunk, my_next_chunk = self._find_and_split_before(before_base1)

        # for each chunk in inserted fragment, add an edge for current fragment
        last_chunk = prev_chunk
        fragment_length = 0
        for chunk in fragment.chunks():
            # also compute how long fragment is
            fragment_length += len(chunk.sequence)
            if last_chunk is None:  # add new chunks at start of fragment
                self.start_chunk = chunk
                self.save()
            else:
                self._assert_not_linked_to(chunk)
                self._add_edges(last_chunk, Edge(from_chunk=last_chunk,
                                                 fragment=self, to_chunk=chunk))
            last_chunk = chunk

        # my_next_chunk may be None, but that's okay, we want to make sure this
        # chunk is the END and not going to be superseded by child fragment
        # appending more chunks!
        self._assert_not_linked_to(my_next_chunk)
        self._add_edges(last_chunk, Edge(from_chunk=last_chunk,
                                         fragment=self, to_chunk=my_next_chunk))

        # shift base_first and base_last for existing chunks
        if before_base1 is not None:
            self.fragment_chunk_location_set.filter(base_first__gte=before_base1)\
                                            .update(base_first=F('base_first')+fragment_length,
                                                    base_last=F('base_last')+fragment_length)

        # add location for new chunks in the new fragment
        values = []
        c = 0
        original_length = self.length
        for chunk in fragment.chunks():
            if before_base1 is not None:
                self.fragment_chunk_location_set.create(
                    chunk=chunk,
                    base_first=before_base1+c,
                    base_last=before_base1+c+len(chunk.sequence)-1
                )
            else:
                self.fragment_chunk_location_set.create(
                    chunk=chunk,
                    base_first=original_length+1+c,
                    base_last=original_length+1+c+len(chunk.sequence)-1
                )
            c += len(chunk.sequence)

    @transaction.atomic()
    def replace_with_fragment(self, before_base1, length_to_remove, fragment):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_fragment(before_base1, fragment)
