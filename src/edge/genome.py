from sqlalchemy.sql import select, and_, func
from edge.models import fragment_table, genome_table, genome_fragment_table
from edge.models import chunk_annotation_table, annotation_table
from edge.models import fragment_chunk_location_table, edge_table, Annotation


class GenomeNotFound(Exception):
    pass


class Genome(object):

    def __init__(self, connector, id):
        self._connector = connector
        self._genome_id = id
        self._genome = self.__get_genome(id)
        self.__updated = []

    def updated(self):
        return self.__updated

    def last_updated(self):
        if len(self.__updated) == 0:
            return None
        return self.__updated[-1]

    def _add_updated(self, operator):
        self.__updated.append(operator)

    @property
    def id(self):
        return self._genome_id

    @property
    def name(self):
        return self._genome['name']

    @property
    def notes(self):
        return self._genome['notes']

    @property
    def parent_id(self):
        return self._genome['parent_id']

    @property
    def parent(self):
        if self._genome['parent_id']:
            return Genome(self._connector, self._genome['parent_id'])
        return None

    @property
    def conn(self):
        return self._connector.conn

    def __get_genome(self, id):
        stmt = \
            select([genome_table.c.name, genome_table.c.parent_id, genome_table.c.notes])\
            .where(genome_table.c.id == id)
        cur = self.conn.execute(stmt)
        genome = None
        for f in cur:
            genome = {'name': f[0], 'parent_id': f[1], 'notes': f[2]}
        if genome is None:
            raise GenomeNotFound('Cannot find genome %s' % (id,))
        return genome

    def fragments(self):
        from edge.fragment import Fragment_Operator

        fragments = []
        stmt = \
            select([fragment_table.c.id])\
            .where(and_(genome_fragment_table.c.genome_id == self.id,
                        genome_fragment_table.c.fragment_id == fragment_table.c.id))
        cur = self.conn.execute(stmt)
        for f in cur:
            fragments.append(Fragment_Operator(self._connector, f[0]))
        return fragments

    def get_fragment_by_id(self, fragment_id):
        f = [x for x in self.fragments() if x.id == fragment_id]
        if len(f) > 0:
            return f[0]
        raise Exception('Fragment with ID %s is not part of the genome' % (fragment_id,))

    def find_annotation(self, name):
        from edge.fragment import Fragment

        stmt = \
            select([annotation_table.c.id,
                    annotation_table.c.name,
                    annotation_table.c.type,
                    annotation_table.c.strand,
                    annotation_table.c.length,
                    chunk_annotation_table.c.annotation_base_first,
                    chunk_annotation_table.c.annotation_base_last,
                    fragment_chunk_location_table.c.base_first,
                    fragment_chunk_location_table.c.base_last,
                    genome_fragment_table.c.fragment_id])\
            .where(and_(
                genome_fragment_table.c.genome_id == self.id,
                genome_fragment_table.c.fragment_id == fragment_chunk_location_table.c.fragment_id,
                chunk_annotation_table.c.chunk_id == fragment_chunk_location_table.c.chunk_id,
                chunk_annotation_table.c.annotation_id == annotation_table.c.id,
                annotation_table.c.name == name))\
            .order_by(annotation_table.c.id,
                      genome_fragment_table.c.fragment_id,
                      fragment_chunk_location_table.c.base_first)

        cur = self.conn.execute(stmt)
        fragments = {}

        for f in cur:
            fragment_id = f[9]
            if fragment_id not in fragments:
                fragments[fragment_id] = []
            annotations = fragments[fragment_id]

            if len(annotations) > 0 and\
               annotations[-1].annotation_id == f[0] and\
               annotations[-1].annotation_last_bp == f[5]-1 and\
               annotations[-1].last_bp == f[7]-1:
                # merge annotation
                annotations[-1].last_bp = f[8]
                annotations[-1].annotation_last_bp = f[6]
            else:
                annotations.append(Annotation(first_bp=f[7], last_bp=f[8],
                                              annotation_id=f[0],
                                              name=f[1], type=f[2], strand=f[3],
                                              annotation_full_length=f[4],
                                              annotation_first_bp=f[5],
                                              annotation_last_bp=f[6]))
        return fragments

    def _edit(self, name=None, notes=None):
        from edge.genome_updater import Genome_Updater
        new_connector = self._connector.new_connector(with_transaction=True)
        u = Genome_Updater(None, new_connector, self.id)
        if notes is not None:
            u._genome['notes'] = notes
        if name is not None:
            u._genome['name'] = name
        return u

    def annotate(self, name=None, notes=None):
        return self._edit(name, notes)

    def update(self, name=None, notes=None):
        from edge.genome_updater import Genome_Updater
        new_connector = self._connector.new_connector(with_transaction=True)
        new_genome_id = Genome_Updater.add_genome(new_connector.conn, self.name, self)
        u = Genome_Updater(self, new_connector, new_genome_id)
        if notes is not None:
            u._genome['notes'] = notes
        if name is not None:
            u._genome['name'] = name
        return u

    def changes(self):
        if self.parent is None:
            return []

        from edge.fragment import Fragment, Fragment_Chunk

        stmt = \
            select([edge_table.c.from_chunk_id,
                    edge_table.c.fragment_id,
                    edge_table.c.to_chunk_id])\
            .where(and_(
                genome_fragment_table.c.genome_id == self.id,
                edge_table.c.fragment_id == genome_fragment_table.c.fragment_id,
                genome_fragment_table.c.inherited == 0,
                fragment_chunk_location_table.c.fragment_id == edge_table.c.fragment_id,
                fragment_chunk_location_table.c.chunk_id == edge_table.c.from_chunk_id,))\
            .order_by(fragment_chunk_location_table.c.fragment_id,
                      fragment_chunk_location_table.c.base_first)

        cur = self.conn.execute(stmt)

        chunks = []
        chunk_ids = []
        for r in cur:
            f = Fragment(self._connector, int(r[1]))
            if r[0] and r[0] not in chunk_ids:
                chunks.append(f.get_chunk(int(r[0])))
                chunk_ids.append(r[0])
            if r[2] and r[2] not in chunk_ids:
                chunks.append(f.get_chunk(int(r[2])))
                chunk_ids.append(r[2])

        return chunks

    def changed_locations_by_fragment(self):
        changes = {}

        for c in self.changes():
            if c.fragment.id not in changes:
                changes[c.fragment.id] = []
            v = changes[c.fragment.id]
            if len(v) == 0 or v[-1][1]+1 != c.location[0]:
                v.append([c.location[0], c.location[1]])
            else:
                v[-1][1] = c.location[1]

        changes = {self.get_fragment_by_id(f): v for f, v in changes.iteritems()}
        return changes
