from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, BigInteger, Integer, String, Index, Boolean, UnicodeText
from sqlalchemy.dialects import sqlite
from sqlalchemy.types import BLOB
from sqlalchemy.schema import Table

Base = declarative_base()

BigIntOrInt = BigInteger()
BigIntOrInt = BigIntOrInt.with_variant(sqlite.INTEGER(), 'sqlite')

chunk_table = Table(
    'chunks', Base.metadata,
    Column('id', BigIntOrInt, primary_key=True, autoincrement=True),
    Column('initial_fragment_id', Integer, nullable=False),
    Column('sequence', String, nullable=True),
    # stores out edges, each edge is (edge ID, fragment ID)
)

fragment_table = Table(
    'fragments', Base.metadata,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('circular', Boolean, nullable=False),
    Column('name', String, nullable=False),
    Column('parent_id', Integer, nullable=True),
    Column('start_chunk_id', BigIntOrInt, nullable=True),
)

edge_table = Table(
    'edges', Base.metadata,
    Column('from_chunk_id', BigIntOrInt, nullable=False),
    Column('fragment_id', Integer, nullable=False),
    # can be null, so we can supersede an edge from a child fragment
    Column('to_chunk_id', BigIntOrInt, nullable=True),
)

annotation_table = Table(
    'annotations', Base.metadata,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('name', String, nullable=False),
    Column('type', String, nullable=False),
    Column('strand', Integer, nullable=True),
    Column('length', Integer, nullable=False),
)

chunk_annotation_table = Table(
    'chunk_annotations', Base.metadata,
    Column('chunk_id', BigIntOrInt, index=True),
    Column('annotation_id', Integer, nullable=False),
    Column('annotation_base_first', Integer, nullable=False),
    Column('annotation_base_last', Integer, nullable=False),
)

genome_table = Table(
    'genomes', Base.metadata,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('name', String, nullable=False),
    Column('parent_id', Integer, nullable=True),
    Column('notes', UnicodeText, nullable=True),
)

genome_fragment_table = Table(
    'genome_fragments', Base.metadata,
    Column('genome_id', Integer, index=True, nullable=False),
    Column('fragment_id', Integer, nullable=False),
    Column('inherited', Boolean, nullable=False),
)

fragment_chunk_location_table = Table(
    'fragment_chunk_locations', Base.metadata,
    Column('fragment_id', Integer, nullable=False),
    Column('chunk_id', BigIntOrInt, nullable=False, index=True),
    Column('base_first', Integer, nullable=False),
    Column('base_last', Integer, nullable=False),
)

Index('idx_fcbf',
      fragment_chunk_location_table.c.fragment_id,
      fragment_chunk_location_table.c.chunk_id,
      fragment_chunk_location_table.c.base_first,
      unique=True)

Index('idx_fcbl',
      fragment_chunk_location_table.c.fragment_id,
      fragment_chunk_location_table.c.chunk_id,
      fragment_chunk_location_table.c.base_last,
      unique=True)


class Edge(object):

    def __init__(self, from_chunk_id, fragment_id, to_chunk_id):
        self.from_chunk_id = from_chunk_id
        self.fragment_id = fragment_id
        self.to_chunk_id = to_chunk_id


class Annotation(object):

    def __init__(self, first_bp, last_bp, annotation_id, name, type, strand,
                 annotation_full_length, annotation_first_bp, annotation_last_bp):
        self.first_bp = first_bp
        self.last_bp = last_bp
        self.annotation_id = annotation_id
        self.name = name
        self.type = type
        self.strand = strand
        self.annotation_full_length = annotation_full_length
        self.annotation_first_bp = annotation_first_bp
        self.annotation_last_bp = annotation_last_bp

    def __str__(self):
        s = []
        if self.annotation_first_bp != 1 or self.annotation_last_bp != self.annotation_full_length:
            s.append('%s (%s-%s)' % (self.name, self.annotation_first_bp, self.annotation_last_bp))
        else:
            s.append(self.name)
        s.append(self.type)
        if self.strand == 1:
            s.append('+')
        else:
            s.append('-')
        return ', '.join(s)
