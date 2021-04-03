from django.db.backends.signals import connection_created

__version__ = '2.17.0'


import random


def get_random_sequence(n):
    return "".join(["agct"[random.randint(0, 3)] for i in range(n)])


def import_gff(name, fn):
    """
    Creates a new genome using the specified GFF file.

    name: Name of genome
    fn: path to GFF file
    """

    from edge.models import Genome
    if Genome.objects.filter(name=name).count() > 0:
        raise Exception('There is already a genome named "%s"' % (name,))
    g = Genome.import_gff(name, fn)
    return g


def _setup_sqlite3(sender, connection, **kwargs):
    """
    Setup SQLite to allow FKs and use in-memory journals
    """

    if connection.vendor == 'sqlite':
        cursor = connection.cursor()
        cursor.execute('PRAGMA foreign_keys = ON;')
        cursor.execute('PRAGMA journal_mode = MEMORY;')


connection_created.connect(_setup_sqlite3)
