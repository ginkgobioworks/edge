from django.db.models import F
from edge.models.chunk import (
    Chunk,
    Chunk_Feature,
    Edge,
    Fragment_Chunk_Location,
)


class Fragment_Writer(object):
    """
    Mixin that includes helpers for updating a fragment.
    """

    def _create_chunk_annotation(self, chunk, feature, feature_base_first, feature_base_last):
        return Chunk_Feature(
            chunk=chunk,
            feature=feature,
            feature_base_first=feature_base_first,
            feature_base_last=feature_base_last,
        )

    def _annotate_chunk(self, chunk, feature, feature_base_first, feature_base_last):
        self._create_chunk_annotation(
            chunk, feature, feature_base_first, feature_base_last
        ).save()

    def _add_reference_chunk(self, start, end, fragment):
        c = Chunk(
            initial_fragment=fragment,
            ref_start_index=start,
            ref_end_index=end
        )
        c.save()
        return c

    def _add_chunk(self, sequence, fragment):
        c = Chunk(sequence=sequence, initial_fragment=fragment)
        c.save()
        return c

    def __reset_chunk_reference(self, chunk, start, end):
        chunk.ref_start_index = start
        chunk.ref_end_index = end
        chunk.save()
        Chunk_Feature.objects.filter(chunk=chunk).delete()
        # if chunk is the start chunk, update the object reference
        if self.start_chunk.id == chunk.id:
            self.start_chunk = chunk
            self.save()

    def __reset_chunk_sequence(self, chunk, sequence):
        chunk.sequence = sequence
        chunk.save()
        Chunk_Feature.objects.filter(chunk=chunk).delete()
        # if chunk is the start chunk, update the object reference
        if self.start_chunk.id == chunk.id:
            self.start_chunk = chunk
            self.save()

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
        cfs = []
        for a in annotations:
            if a.feature.strand > 0:
                a1 = (
                    a.feature,
                    a.feature_base_first,
                    a.feature_base_first + bps_to_split - 1,
                )
                a2 = (
                    a.feature,
                    a.feature_base_first + bps_to_split,
                    a.feature_base_last
                )
            else:
                a1 = (
                    a.feature,
                    a.feature_base_last - bps_to_split + 1,
                    a.feature_base_last
                )
                a2 = (
                    a.feature,
                    a.feature_base_first,
                    a.feature_base_last - bps_to_split
                )
            cfs.append(self._create_chunk_annotation(split1, *a1))
            cfs.append(self._create_chunk_annotation(split2, *a2))
        Chunk_Feature.bulk_create(cfs)

    # make sure you call this atomically! otherwise we may have corrupted chunk
    # and index
    def __split_chunk(self, chunk, bps_to_split):
        # retrieve relevant fragments
        candidates = (
            chunk.fragment_chunk_location_set.filter(
                fragment__fragment_index__fresh=True
            )
            .values_list("fragment_id")
            .distinct()
        )
        index_to_update = [row[0] for row in candidates]

        # splitted chunk should be "created" by the fragment that created the
        # original chunk
        if chunk.is_sequence_based:
            chunk_sequence = chunk.get_sequence()
            s1 = chunk_sequence[0:bps_to_split]
            s2 = chunk_sequence[bps_to_split:]
            split2 = self._add_chunk(s2, chunk.initial_fragment)
            self.__reset_chunk_sequence(chunk, s1)
        elif chunk.is_reference_based:
            # note reference based chunks have 1-based indexing
            split_start = chunk.ref_start_index + bps_to_split
            split2 = self._add_reference_chunk(
                split_start, chunk.ref_end_index, chunk.initial_fragment
            )
            self.__reset_chunk_reference(chunk, chunk.ref_start_index, split_start - 1)

        # move all edges from original chunk to second chunk
        Edge.objects.filter(from_chunk=chunk).update(from_chunk=split2)
        # verify chunk edges been updated
        assert chunk.out_edges.count() == 0

        # add edge from chunk to split2
        self._add_edges(
            chunk,
            Edge(from_chunk=chunk, fragment=chunk.initial_fragment, to_chunk=split2),
        )

        # add location for new chunk - note that the following code has to
        # update the location index for each of the fragment that uses the
        # chunk. there may be a lot of fragments, potentially a big scalability
        # problem.
        entries = []
        for fcl in chunk.fragment_chunk_location_set.filter(
            fragment_id__in=index_to_update
        ):
            entries.append(
                Fragment_Chunk_Location(
                    fragment_id=fcl.fragment_id,
                    chunk_id=split2.id,
                    base_first=fcl.base_first + bps_to_split,
                    base_last=fcl.base_last,
                )
            )
        Fragment_Chunk_Location.bulk_create(entries)

        # adjust chunk location index for existing chunk
        chunk.fragment_chunk_location_set.filter(
            fragment_id__in=index_to_update
        ).update(base_last=F("base_first") + bps_to_split - 1)
        return split2

    def _find_chunk_prev_next(self, before_base1):
        prev_chunk = None
        next_chunk = None
        chunk = None
        bases_visited = None

        if before_base1 is not None:
            q = self.fragment_chunk_location_set.filter(
                base_first__lte=before_base1, base_last__gte=before_base1
            )
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
                prev_chunk = self.fragment_chunk_location_set.get(
                    base_last=total_bases
                ).chunk
            bases_visited = total_bases

        return prev_chunk, chunk, next_chunk, bases_visited

    def _find_and_split_before(self, before_base1):

        if before_base1 is not None and before_base1 <= 0:
            raise Exception("chunk index should be 1-based")

        chunk = None
        prev_chunk = None
        next_chunk = None
        bases_visited = 0

        prev_chunk, chunk, next_chunk, bases_visited = self._find_chunk_prev_next(
            before_base1
        )

        # found the bp we are looking for
        if before_base1 is not None and bases_visited >= before_base1:
            chunk_len = chunk.length

            # can avoid splitting if first bp in this chunk is before_base1
            if bases_visited - chunk_len + 1 == before_base1:
                # print('no need to split before %s' % (before_base1,))
                return prev_chunk, chunk

            # otherwise, have to split the chunk
            first_bp_in_chunk = bases_visited - chunk_len + 1
            bps_to_split = before_base1 - first_bp_in_chunk

            # save original annotations, which will be trashed in __split_chunk
            # (which calls __reset_chunk_sequence)
            original_annotations = self.fragment_chunk(chunk).annotations()
            # split chunk
            split2 = self.__split_chunk(chunk, bps_to_split)
            chunk = chunk.reload()

            # split up annotations as well
            self._split_annotations(original_annotations, bps_to_split, chunk, split2)
            return chunk, split2

        else:  # found end of the fragment
            return prev_chunk, None
