from sqlalchemy.sql import select, and_
from edge.models import fragment_table, chunk_table, annotation_table, chunk_annotation_table
from edge.models import fragment_chunk_location_table, edge_table, Edge
from edge.fragment import Fragment, Fragment_Operator


class Fragment_Writer(Fragment):

    def __init__(self, operator, *args, **kwargs):
        super(Fragment_Writer, self).__init__(*args, **kwargs)
        self.__parent_operator = operator

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if type is None:  # only commit when there are no exceptions
            f = self.commit()
            if self.__parent_operator:
                self.__parent_operator._add_updated(f)

    def can_use_location_index(self):
        return False

    def commit(self):
        stmt = \
            fragment_table.update()\
            .where(fragment_table.c.id == self._fragment_id)\
            .values(start_chunk_id=self.start_chunk_id)
        self.conn.execute(stmt)
        self._connector.commit()
        return Fragment_Operator(self._connector, self._fragment_id)

    def _add_annotation(self, name, type, length, strand):
        if strand not in (1, -1, None):
            raise Exception('Strand must be 1, -1, or None')
        stmt = annotation_table.insert().values(name=name, type=type, length=length, strand=strand)
        res = self.conn.execute(stmt)
        return res.inserted_primary_key[0]

    def _annotate_chunk(self, chunk_id, annotation_id, annotation_base_first, annotation_base_last):
        stmt = chunk_annotation_table.insert().values(
            chunk_id=chunk_id,
            annotation_id=annotation_id,
            annotation_base_first=annotation_base_first,
            annotation_base_last=annotation_base_last
        )
        res = self.conn.execute(stmt)

    def _add_chunk(self, sequence, fragment_id):
        stmt = chunk_table.insert().values(sequence=sequence, initial_fragment_id=fragment_id)
        res = self.conn.execute(stmt)
        return res.inserted_primary_key[0]

    def _reset_chunk_sequence(self, chunk_id, sequence):
        stmt = chunk_table.update().where(chunk_table.c.id == chunk_id).values(sequence=sequence)
        res = self.conn.execute(stmt)
        stmt = chunk_annotation_table.delete().where(chunk_annotation_table.c.chunk_id == chunk_id)
        res = self.conn.execute(stmt)

    def _assert_not_linked_to(self, chunk_id):
        stmt = \
            select([edge_table.c.from_chunk_id])\
            .where(and_(edge_table.c.to_chunk_id == chunk_id,
                        edge_table.c.fragment_id == self.id))
        res = self.conn.execute(stmt)
        res = [r for r in res]
        if len(res) > 0:
            raise Exception('Fragment %s already linked to chunk %s' % (self.id, chunk_id))

    def _add_edges(self, chunk_id, *edges):
        # fetch previous list of out edges
        stmt = \
            select([edge_table.c.from_chunk_id,
                    edge_table.c.fragment_id,
                    edge_table.c.to_chunk_id])\
            .where(edge_table.c.from_chunk_id == chunk_id)
        res = self.conn.execute(stmt)
        existing_edges = [Edge(**r) for r in res]
        res.close()

        # get list of fragment IDs in new edges
        new_fragment_ids = [e.fragment_id for e in edges]

        # remove old edges with same fragment ids as new edges
        for e in existing_edges:
            if e.fragment_id in new_fragment_ids:
                stmt = \
                    edge_table.delete()\
                    .where(and_(edge_table.c.from_chunk_id == e.from_chunk_id,
                                edge_table.c.fragment_id == e.fragment_id,
                                edge_table.c.to_chunk_id == e.to_chunk_id))
                self.conn.execute(stmt)

        # add new edges
        stmt = edge_table.insert()
        values = []
        for edge in edges:
            values.append(dict(from_chunk_id=edge.from_chunk_id,
                               fragment_id=edge.fragment_id,
                               to_chunk_id=edge.to_chunk_id))
        self.conn.execute(stmt, values)

    def _split_annotations(self, annotations, bps_to_split, split1, split2):
        for annotation in annotations:
            a1 = (annotation.annotation_id,
                  annotation.annotation_first_bp,
                  annotation.annotation_first_bp+bps_to_split-1)
            a2 = (annotation.annotation_id,
                  annotation.annotation_first_bp+bps_to_split,
                  annotation.annotation_last_bp)
            self._annotate_chunk(split1, *a1)
            self._annotate_chunk(split2, *a2)

    def _split_chunk(self, chunk, s1, s2):
        # splitted chunk should be "created" by the fragment that created the
        # original chunk
        split2 = self._add_chunk(s2, chunk.initial_fragment_id)
        self._reset_chunk_sequence(chunk.id, s1)

        # move all edges from original chunk to second chunk
        stmt = \
            edge_table.update()\
            .where(edge_table.c.from_chunk_id == chunk.id).values(from_chunk_id=split2)
        res = self.conn.execute(stmt)
        chunk.reload()

        # add edge from chunk to split2
        self._add_edges(chunk.id, Edge(chunk.id, chunk.initial_fragment_id, split2))

        # add location for new chunk
        stmt = \
            select([fragment_chunk_location_table.c.fragment_id,
                    fragment_chunk_location_table.c.chunk_id,
                    fragment_chunk_location_table.c.base_first,
                    fragment_chunk_location_table.c.base_last])\
            .where(fragment_chunk_location_table.c.chunk_id == chunk.id)

        for r in self.conn.execute(stmt):
            fragment_id = r[0]
            chunk_id = r[1]
            bp_first = r[2]
            bp_last = r[3]

            insert_stmt = \
                fragment_chunk_location_table.insert()\
                .values(fragment_id=fragment_id,
                        chunk_id=split2,
                        base_first=bp_first+len(s1),
                        base_last=bp_last)
            self.conn.execute(insert_stmt)

        # adjust chunk location index for existing chunk
        update_stmt = fragment_chunk_location_table.update()\
            .where(fragment_chunk_location_table.c.chunk_id == chunk.id)\
            .values(base_last=fragment_chunk_location_table.c.base_first+len(s1)-1)
        self.conn.execute(update_stmt)

        return split2

    def _find_chunk_prev_next_by_walking(self, before_base1):
        prev_chunk_id = None
        next_chunk_id = None
        chunk = None
        bases_visited = 0

        chunk_id = self.start_chunk_id
        while chunk_id is not None:
            chunk = self.get_chunk(chunk_id)
            sequence = chunk.sequence
            next_chunk_id = chunk.next_chunk_id
            bases_visited += len(chunk.sequence)
            if before_base1 is not None and bases_visited >= before_base1:
                break
            prev_chunk_id = chunk_id
            chunk_id = next_chunk_id
            next_chunk_id = None
            chunk = None

        return prev_chunk_id, chunk, next_chunk_id, bases_visited

    def _find_chunk_prev_next_by_location_index(self, before_base1):
        prev_chunk_id = None
        next_chunk_id = None
        chunk = None
        bases_visited = None

        if before_base1 is not None:
            stmt = \
                select([fragment_chunk_location_table.c.chunk_id,
                        fragment_chunk_location_table.c.base_first,
                        fragment_chunk_location_table.c.base_last])\
                .where(and_(fragment_chunk_location_table.c.fragment_id == self.id,
                            fragment_chunk_location_table.c.base_first <= before_base1,
                            fragment_chunk_location_table.c.base_last >= before_base1))

            res = self.conn.execute(stmt)
            res = [r for r in res]
            if len(res) > 0:
                chunk_id = res[0][0]
                chunk_base_first = res[0][1]
                chunk_base_last = res[0][2]
                chunk = self.get_chunk(chunk_id)
                bases_visited = chunk_base_last
                next_chunk_id = chunk.next_chunk_id

                # find prev chunk
                if chunk_base_first == 1:
                    prev_chunk_id = None
                else:
                    stmt = \
                        select([fragment_chunk_location_table.c.chunk_id])\
                        .where(and_(
                            fragment_chunk_location_table.c.fragment_id == self.id,
                            fragment_chunk_location_table.c.base_last == chunk_base_first-1
                        ))
                    res = self.conn.execute(stmt)
                    prev_chunk_id = res.fetchone()[0]

        if chunk is None:  # after all sequence ended, need last chunk and total bases
            total_bases = self.length
            stmt = \
                select([fragment_chunk_location_table.c.chunk_id])\
                .where(and_(fragment_chunk_location_table.c.fragment_id == self.id,
                            fragment_chunk_location_table.c.base_last == total_bases))
            res = self.conn.execute(stmt)
            prev_chunk_id = res.fetchone()[0]
            bases_visited = total_bases

        return prev_chunk_id, chunk, next_chunk_id, bases_visited

    def _find_and_split_before(self, before_base1):

        if before_base1 is not None and before_base1 <= 0:
            raise Exception('chunk index should be 1-based')

        chunk_id = None
        chunk = None
        bases_visited = 0
        prev_chunk_id = None
        next_chunk_id = None
        sequence = None

        if self.can_use_location_index():
            prev_chunk_id, chunk, next_chunk_id, bases_visited = \
                self._find_chunk_prev_next_by_location_index(before_base1)
            if chunk:
                sequence = chunk.sequence
                chunk_id = chunk.id

        else:
            prev_chunk_id, chunk, next_chunk_id, bases_visited = \
                self._find_chunk_prev_next_by_walking(before_base1)
            if chunk:
                sequence = chunk.sequence
                chunk_id = chunk.id

        # found the bp we are looking for
        if before_base1 is not None and bases_visited >= before_base1:

            # can avoid splitting if first bp in this chunk is before_base1
            if bases_visited-len(sequence)+1 == before_base1:
                #print 'no need to split before %s, which contains %s' % (before_base1, sequence)
                return prev_chunk_id, chunk_id

            # otherwise, have to split the chunk
            first_bp_in_chunk = bases_visited-len(sequence)+1
            bps_to_split = before_base1-first_bp_in_chunk
            s1 = sequence[0:bps_to_split]
            s2 = sequence[bps_to_split:]

            # save original annotations, which will be trashed in _split_chunk (which
            # calls _reset_chunk_sequence)
            original_annotations = chunk.annotations()
            # split chunk
            split2 = self._split_chunk(chunk, s1, s2)
            # split up annotations as well
            self._split_annotations(original_annotations, bps_to_split, chunk_id, split2)
            return chunk_id, split2

        else:  # found end of the fragment
            return prev_chunk_id, None


class Fragment_Annotator(Fragment_Writer):

    def can_use_location_index(self):
        return True

    def annotate(self, first_base1, last_base1, annotation_name, annotation_type, strand):

        if self.circular and last_base1 < first_base1:
            # has to figure out the total length from last chunk
            length = self.length-first_base1+1+last_base1
        else:
            length = last_base1-first_base1+1
            if length <= 0:
                raise Exception('Annotation must have length one or more')

        prev_chunk_id, annotation_start = self._find_and_split_before(first_base1)
        annotation_end, next_chunk_id = self._find_and_split_before(last_base1+1)

        new_annotation_id = self._add_annotation(annotation_name, annotation_type, length, strand)

        # now, starting with chunk annotation_start, walk through chunks until we
        # hit annotation_end, and add annotation for each chunk
        chunk = self.get_chunk(annotation_start)
        a_i = 1
        while True:
            self._annotate_chunk(chunk.id, new_annotation_id, a_i, a_i+len(chunk.sequence)-1)
            a_i += len(chunk.sequence)
            if chunk.id == annotation_end:
                break
            if chunk.next_chunk_id:
                chunk = self.get_chunk(chunk.next_chunk_id)
            else:
                chunk = self.get_chunk(self.start_chunk_id)


class Fragment_Updater(Fragment_Writer):

    @staticmethod
    def add_fragment(conn, name, circular, parent_id=None, start_chunk_id=None):
        stmt = fragment_table.insert().values(
            name=name, circular=circular, parent_id=parent_id, start_chunk_id=start_chunk_id
        )
        res = conn.execute(stmt)
        return res.inserted_primary_key[0]

    def __init__(self, operator, use_location_index, *args, **kwargs):
        super(Fragment_Updater, self).__init__(operator, *args, **kwargs)
        self.__use_location_index = use_location_index

    def can_use_location_index(self):
        return self.__use_location_index

    def _index_fragment_chunk_locations(self):
        # remove old index
        stmt = \
            fragment_chunk_location_table.delete()\
            .where(fragment_chunk_location_table.c.fragment_id == self.id)
        self.conn.execute(stmt)

        # go through each chunk and add new index
        i = 1
        for chunk in self.chunks():
            if len(chunk.sequence) > 0:
                stmt = fragment_chunk_location_table.insert().values(
                    fragment_id=self.id, chunk_id=chunk.id,
                    base_first=i, base_last=i+len(chunk.sequence)-1,
                )
                self.conn.execute(stmt)
                i += len(chunk.sequence)

    def commit(self):
        # add fragment chunk location index if we have not been keeping it
        if not self.can_use_location_index():
            self._index_fragment_chunk_locations()
        return super(Fragment_Updater, self).commit()

    def insert_bases(self, before_base1, sequence):
        # find chunks before and containing the insertion point
        prev_chunk_id, chunk_id = self._find_and_split_before(before_base1)

        # create new chunk
        new_chunk = self._add_chunk(sequence, self._fragment_id)

        if prev_chunk_id is not None:  # add chunks after prev_chunk_id
            self._add_edges(prev_chunk_id, Edge(prev_chunk_id, self._fragment_id, new_chunk))

        else:  # add chunks at start of fragment
            self._set_start_chunk_id(new_chunk)

        # chunk_id may be None, but that's okay, we want to make sure this chunk is
        # the END and not going to be superseded by child fragment appending more
        # chunks!
        self._add_edges(new_chunk, Edge(new_chunk, self._fragment_id, chunk_id))

        # update location index
        if self.can_use_location_index():

            # shift base_first and base_last for existing chunks
            if before_base1 is not None:
                update_stmt = \
                    fragment_chunk_location_table.update()\
                    .where(and_(fragment_chunk_location_table.c.fragment_id == self._fragment_id,
                                fragment_chunk_location_table.c.base_first >= before_base1))\
                    .values(base_first=fragment_chunk_location_table.c.base_first+len(sequence),
                            base_last=fragment_chunk_location_table.c.base_last+len(sequence))
                self.conn.execute(update_stmt)

            # insert location for new chunk
            if before_base1 is not None:
                insert_stmt = \
                    fragment_chunk_location_table.insert()\
                    .values(fragment_id=self._fragment_id, chunk_id=new_chunk,
                            base_first=before_base1, base_last=before_base1+len(sequence)-1)
                self.conn.execute(insert_stmt)
            else:
                fragment_length = self.length
                insert_stmt = \
                    fragment_chunk_location_table.insert()\
                    .values(fragment_id=self._fragment_id, chunk_id=new_chunk,
                            base_first=fragment_length+1,
                            base_last=fragment_length+1+len(sequence)-1)
                self.conn.execute(insert_stmt)

    def remove_bases(self, before_base1, length):
        if length <= 0:
            raise Exception('Cannot remove less than one base pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        prev_chunk_id, removal_start = self._find_and_split_before(before_base1)
        removal_end, next_chunk_id = self._find_and_split_before(before_base1+length)

        if prev_chunk_id is not None:  # remove chunks after prev_chunk_id
            if next_chunk_id:  # delete prior to end of fragment
                self._add_edges(prev_chunk_id,
                                Edge(prev_chunk_id, self._fragment_id, next_chunk_id))
            else:
                # delete all remaining chunks in fragment by adding an edge with None
                # as target, superseding any child that may have added chunks after
                # prev_chunk_id
                self._add_edges(prev_chunk_id,
                                Edge(prev_chunk_id, self._fragment_id, None))

        else:  # remove chunks at start of fragment
            if next_chunk_id is None:
                raise Exception('Cannot remove entire fragment')
            self._set_start_chunk_id(next_chunk_id)

        # update location index
        if self.can_use_location_index():

            # remove location for deleted chunks
            remove_stmt = \
                fragment_chunk_location_table.delete()\
                .where(and_(fragment_chunk_location_table.c.fragment_id == self._fragment_id,
                            fragment_chunk_location_table.c.base_first >= before_base1,
                            fragment_chunk_location_table.c.base_first < before_base1+length))
            self.conn.execute(remove_stmt)

            # shift base_first and base_last for existing chunks
            update_stmt = \
                fragment_chunk_location_table.update()\
                .where(and_(fragment_chunk_location_table.c.fragment_id == self._fragment_id,
                            fragment_chunk_location_table.c.base_first > before_base1))\
                .values(base_first=fragment_chunk_location_table.c.base_first-length,
                        base_last=fragment_chunk_location_table.c.base_last-length)
            self.conn.execute(update_stmt)

    def replace_bases(self, before_base1, length_to_remove, sequence):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_bases(before_base1, sequence)

    def insert_fragment(self, before_base1, fragment):

        # find chunks before and containing the insertion point
        prev_chunk_id, my_next_chunk_id = self._find_and_split_before(before_base1)

        # for each chunk in inserted fragment, add an edge for current fragment
        last_chunk_id = prev_chunk_id
        fragment_length = 0
        for chunk in fragment.chunks():
            # also compute how long fragment is
            fragment_length += len(chunk.sequence)
            if last_chunk_id is None:  # add new chunks at start of fragment
                self._set_start_chunk_id(chunk.id)
            else:
                self._assert_not_linked_to(chunk.id)
                self._add_edges(last_chunk_id, Edge(last_chunk_id, self._fragment_id, chunk.id))
            last_chunk_id = chunk.id

        # my_next_chunk_id may be None, but that's okay, we want to make sure this
        # chunk is the END and not going to be superseded by child fragment
        # appending more chunks!
        self._assert_not_linked_to(my_next_chunk_id)
        self._add_edges(last_chunk_id, Edge(last_chunk_id, self._fragment_id, my_next_chunk_id))

        # update location index
        if self.can_use_location_index():

            # shift base_first and base_last for existing chunks
            if before_base1 is not None:
                update_stmt = fragment_chunk_location_table.update()\
                    .where(and_(fragment_chunk_location_table.c.fragment_id == self._fragment_id,
                                fragment_chunk_location_table.c.base_first >= before_base1))\
                    .values(base_first=fragment_chunk_location_table.c.base_first+fragment_length,
                            base_last=fragment_chunk_location_table.c.base_last+fragment_length)
                self.conn.execute(update_stmt)

            # add location for new chunks in the new fragment
            values = []
            c = 0
            for chunk in fragment.chunks():
                if before_base1 is not None:
                    values.append(dict(fragment_id=self._fragment_id,
                                       chunk_id=chunk.id,
                                       base_first=before_base1+c,
                                       base_last=before_base1+c+len(chunk.sequence)-1))
                else:
                    values.append(dict(fragment_id=self._fragment_id,
                                       chunk_id=chunk.id,
                                       base_first=self.length+1+c,
                                       base_last=self.length+1+c+len(chunk.sequence)-1))
                c += len(chunk.sequence)
            insert_stmt = fragment_chunk_location_table.insert()
            self.conn.execute(insert_stmt, values)

    def replace_with_fragment(self, before_base1, length_to_remove, fragment):

        if length_to_remove <= 0:
            raise Exception('Cannot remove less than one chunk pair')
        if before_base1 is None:
            raise Exception('Missing position to remove sequences')

        self.remove_bases(before_base1, length_to_remove)
        self.insert_fragment(before_base1, fragment)
