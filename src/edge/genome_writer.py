from contextlib import contextmanager
from sqlalchemy.sql import and_
from BCBio import GFF
from edge.models import genome_fragment_table
from edge.connector import Shared_Connector
from edge.models import genome_table
from edge.genome import Genome


class Genome_Updater(Genome):

    @staticmethod
    def add_genome(conn, name, parent=None, notes=None):
        if parent:
            stmt = genome_table.insert().values(name=name, notes=notes, parent_id=parent.id)
        else:
            stmt = genome_table.insert().values(name=name, notes=notes, parent_id=None)
        res = conn.execute(stmt)
        new_genome_id = res.inserted_primary_key[0]
        if parent:
            stmt = genome_fragment_table.insert()
            values = []
            for f in parent.fragments():
                values.append(dict(genome_id=new_genome_id, fragment_id=f.id, inherited=True))
            if len(values) > 0:
                conn.execute(stmt, values)
        return new_genome_id

    def __init__(self, genome, *args, **kwargs):
        super(Genome_Updater, self).__init__(*args, **kwargs)
        self.__parent_genome = genome

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if type is None:  # only commit when there are no exceptions
            f = self.commit()
            if self.__parent_genome:
                self.__parent_genome._add_updated(f)

    def commit(self):
        stmt = \
            genome_table.update()\
            .where(genome_table.c.id == self.id)\
            .values(notes=self.notes, name=self.name)
        self.conn.execute(stmt)
        self._connector.commit()
        return Genome(self._connector, self.id)

    @contextmanager
    def annotate_fragment_by_name(self, name):
        f = [x for x in self.fragments() if x.name == name]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        u = f[0].annotate()
        yield u
        u.commit()

    @contextmanager
    def annotate_fragment_by_fragment_id(self, fragment_id):
        f = [x for x in self.fragments() if x.id == fragment_id]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        u = f[0].annotate()
        yield u
        u.commit()

    @contextmanager
    def update_fragment_by_name(self, name):
        if self.__parent_genome is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments() if x.name == name]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        u = f[0].update(name)
        yield u
        self._add_updated_fragment(u.commit())

    @contextmanager
    def update_fragment_by_fragment_id(self, fragment_id):
        if self.__parent_genome is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments() if x.id == fragment_id]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        u = f[0].update(f[0].name)
        yield u
        self._add_updated_fragment(u.commit())

    def add_fragment(self, name, sequence, circular=False, annotate=True):
        from edge.fragment_updater import Fragment_Updater, Fragment_Annotator

        if len(sequence) == 0:
            raise Exception('Cannot create a fragment of length zero')

        fragment_id = Fragment_Updater.add_fragment(self.conn, name, circular, None, None)
        u = Fragment_Updater(None, False,
                             self._connector.new_connector(with_transaction=True), fragment_id)
        u.insert_bases(None, sequence)
        f = u.commit()

        stmt = \
            genome_fragment_table.insert()\
            .values(genome_id=self.id, fragment_id=fragment_id, inherited=False)
        self.conn.execute(stmt)

        if annotate:
            return f.annotate()
        return f

    def _add_updated_fragment(self, fragment):
        existing_fragment_ids = [f.id for f in self.fragments()]
        if fragment.parent_id in existing_fragment_ids:
            stmt = \
                genome_fragment_table.update()\
                .where(and_(genome_fragment_table.c.genome_id == self.id,
                            genome_fragment_table.c.fragment_id == fragment.parent_id))\
                .values(fragment_id=fragment.id, inherited=False)
            self.conn.execute(stmt)
        else:
            raise Exception('Fragment parent not part of the genome')

    def import_gff(self, gff_fasta_fn):
        in_file = gff_fasta_fn
        in_handle = open(in_file)
        for rec in GFF.parse(in_handle):
            f = self.add_fragment(rec.id, str(rec.seq))
            #print '%s: %s' % (rec.id, len(str(rec.seq)))
            for feature in rec.features:
                # skip features that cover the entire sequence
                if feature.location.start == 0 and feature.location.end == len(str(rec.seq)):
                    continue
                name = feature.id
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0]
                elif 'Name' in feature.qualifiers:
                    name = feature.qualifiers['Name'][0]
                # start in Genbank format is start after, so +1 here
                f.annotate(feature.location.start+1,
                           feature.location.end,
                           name,
                           feature.type,
                           feature.strand)
                #print '      %s %s: %s %s' % (feature.type, name, feature.location, feature.strand)
            f.commit()
        in_handle.close()
