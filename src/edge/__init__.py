def import_gff(name, fn):
    """
    Creates a new genome using the specified GFF file.

    name: Name of genome
    fn: path to GFF file
    """

    from edge.models import Genome
    g = Genome.create(name)
    u = g.edit()
    u.import_gff(fn)


def _setup_sqlite3(sender, connection, **kwargs):
    """
    Setup SQLite to allow FKs and use in-memory journals
    """

    if connection.vendor == 'sqlite':
        cursor = connection.cursor()
        cursor.execute('PRAGMA foreign_keys = ON;')
        cursor.execute('PRAGMA journal_mode = MEMORY;')

from django.db.backends.signals import connection_created
connection_created.connect(_setup_sqlite3)
