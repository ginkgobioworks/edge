import os
import tempfile
from edge.models import Base
from sqlalchemy import create_engine


class DB(object):
    """
    Class to manage connection to database.
    """

    @staticmethod
    def __sqlite_path(dbfn):
        """
        Returns sqlite path of DB
        """
        return 'sqlite:///%s' % (dbfn,)

    @staticmethod
    def __create_engine(dbfn):
        """
        Returns sqlalchemy engine to database
        """
        # turn echo=True for logging
        return create_engine(DB.__sqlite_path(dbfn), echo=False)

    def __init__(self, dbfn, engine=None):
        """
        Constructor: takes a db filename and optionally an existing engine.
        """
        if not os.path.exists(dbfn):
            raise Exception('Database file %s does not exist' % (dbfn,))
        self.__path = dbfn
        self.__engine = engine if engine else DB.__create_engine(dbfn)

    @property
    def engine(self):
        """
        Returns engine
        """
        return self.__engine

    @property
    def path(self):
        """
        Returns path of DB
        """
        return self.__path

    @classmethod
    def create(cls, dbfn=None):
        """
        Create a new database. If file is given, will use file as new
        database's file name.
        """

        if dbfn is None:
            with tempfile.NamedTemporaryFile(mode='w+b', delete=False, suffix='.db') as dbf:
                dbf.close()
                dbfn = dbf.name
        try:
            os.unlink(dbfn)
        except OSError:
            pass
        engine = DB.__create_engine(dbfn)
        Base.metadata.create_all(engine)
        return cls(dbfn, engine)
