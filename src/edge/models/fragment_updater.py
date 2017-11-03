from django.db.models import F
from edge.models.chunk import Edge


class Fragment_Updater(object):
    """
    Mixin for sequence manipulation within fragment.
    """

    def _assert_not_linked_to(self, chunk):
        if Edge.objects.filter(to_chunk=chunk, fragment=self).count() > 0:
            raise Exception('Fragment %s already linked to chunk %s' % (self.id, chunk.id))

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
            base_first=cur_fragment_length + 1,
            base_last=cur_fragment_length + 1 + len(sequence) - 1
        )
        return new_chunk

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
                                            .update(base_first=F('base_first') + len(sequence),
                                                    base_last=F('base_last') + len(sequence))

        # insert location for new chunk
        if before_base1 is not None:
            self.fragment_chunk_location_set.create(
                chunk=new_chunk,
                base_first=before_base1,
                base_last=before_base1 + len(sequence) - 1
            )
        else:
            fragment_length = self.length
            self.fragment_chunk_location_set.create(
                chunk=new_chunk,
                base_first=fragment_length + 1,
                base_last=fragment_length + 1 + len(sequence) - 1
            )

    def remove_bases(self, before_base1, length):
        if length <= 0:
            raise Exception('Cannot remove less than one base pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        prev_chunk, dont_use_stale_after_second_split = self._find_and_split_before(before_base1)
        removal_end, next_chunk = self._find_and_split_before(before_base1 + length)

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
                                                base_first__lt=before_base1 + length).delete()

        # shift base_first and base_last for existing chunks
        self.fragment_chunk_location_set.filter(base_first__gt=before_base1)\
                                        .update(base_first=F('base_first') - length,
                                                base_last=F('base_last') - length)

    def replace_bases(self, before_base1, length_to_remove, sequence):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_bases(before_base1, sequence)

    def insert_fragment(self, before_base1, fragment):
        fragment = fragment.indexed_fragment()

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
                                            .update(base_first=F('base_first') + fragment_length,
                                                    base_last=F('base_last') + fragment_length)

        # add location for new chunks in the new fragment
        c = 0
        original_length = self.length
        for chunk in fragment.chunks():
            if before_base1 is not None:
                self.fragment_chunk_location_set.create(
                    chunk=chunk,
                    base_first=before_base1 + c,
                    base_last=before_base1 + c + len(chunk.sequence) - 1
                )
            else:
                self.fragment_chunk_location_set.create(
                    chunk=chunk,
                    base_first=original_length + 1 + c,
                    base_last=original_length + 1 + c + len(chunk.sequence) - 1
                )
            c += len(chunk.sequence)

    def replace_with_fragment(self, before_base1, length_to_remove, fragment):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_fragment(before_base1, fragment)
