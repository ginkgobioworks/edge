from django.db.backends.signals import connection_created


def setup_sqlite3(sender, connection, **kwargs):
    """Enable integrity constraint with sqlite."""
    if connection.vendor == 'sqlite':
        cursor = connection.cursor()
        cursor.execute('PRAGMA foreign_keys = ON;')
        cursor.execute('PRAGMA journal_mode = MEMORY;')

connection_created.connect(setup_sqlite3)
