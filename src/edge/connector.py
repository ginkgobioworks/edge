from sqlalchemy.sql import select
from edge.models import genome_table, fragment_table, genome_fragment_table
from edge.genome import Genome


class Connector(object):
    """
    Represents a connection to a database. If a connector is under
    transactional management, creating a new connector returns a shared
    connection to the database.
    """

    @classmethod
    def create_db(cls, dbfn):
        """
        Create a new database using the specified file name.
        """
        from edge.db import DB
        dbobj = DB.create(dbfn)
        return cls(dbobj.engine)

    @classmethod
    def open_db(cls, dbfn):
        """
        Opens an existing database located at the specified file name.
        """
        from edge.db import DB
        dbobj = DB(dbfn)
        return cls(dbobj.engine)

    def __init__(self, engine, with_transaction=False):
        """
        Creates a new connection, using the sqlalchemy engine provided.
        with_transaction defaults to False. If set to True, starts a
        transaction immediately.
        """

        self.__engine = engine
        self.__c = self.__engine.connect()
        if with_transaction:
            self.__transaction = self.__c.begin()
            self.__committed = False
        else:
            self.__transaction = None
            self.__committed = True

    @property
    def conn(self):
        """
        Returns sqlalchemy connection to database.
        """
        return self.__c

    def close(self):
        """
        Closes connection
        """
        self.__c.close()
        self.__c = None

    def commit(self):
        """
        Commits transaction. Connection is no longer under transaction
        management after this point.
        """
        if self.__transaction:
            #print 'really committing'
            self.__transaction.commit()
            self.__transaction = None
        self.__committed = True
        return self

    def new_connector(self, with_transaction):
        """
        Creates a new connector. Passing with_transaction argument to new
        connector's constructor. If current connection has not been committed,
        then creates a Shared_Connector object instead, which shares the same
        connection as the current connector object. Note that you cannot commit
        from a Shared_Connector object.
        """

        if self.__committed:
            #print 'create new connector'
            return Connector(self.__engine, with_transaction)
        else:
            #print 'sharing connector'
            return Shared_Connector(self)

    def non_genomic_fragments(self):
        """
        Returns list of non genomic fragments.
        """
        from edge.fragment import Fragment_Operator

        # get all fragments
        stmt = select([fragment_table.c.id])
        cur = self.conn.execute(stmt)
        fragment_ids = [row[0] for row in cur]
        # get genomic fragments
        stmt = select([genome_fragment_table.c.fragment_id])
        cur = self.conn.execute(stmt)
        genome_fragment_ids = [row[0] for row in cur]
        # subtract to get non-genomic fragments
        fragments = [Fragment_Operator(self, id)
                     for id in list(set(fragment_ids)-set(genome_fragment_ids))]
        return fragments

    def genomes(self):
        """
        Returns list of genomes from the database.
        """

        res = []
        stmt = select([genome_table.c.id])
        cur = self.conn.execute(stmt)
        for row in cur:
            res.append(Genome(self, row[0]))
        return res

    def create_genome(self, name, notes=None):
        """
        Creates a new genome in the database.
        """

        from edge.genome_updater import Genome_Updater

        new_connector = self.new_connector(with_transaction=True)
        new_genome_id = Genome_Updater.add_genome(new_connector.conn, name, notes=notes)
        new_connector.commit()
        return Genome(new_connector, new_genome_id)

    def create_fragment_with_sequence(self, name, sequence, circular=False):
        """
        Creates a fragment in the database. The fragment will not be associated
        with any genome.
        """

        from edge.fragment_updater import Fragment_Updater

        if len(sequence) == 0:
            raise Exception('Cannot create a fragment of length zero')

        new_connector = self.new_connector(with_transaction=True)
        fragment_id = Fragment_Updater.add_fragment(new_connector.conn, name, circular, None, None)
        updater = Fragment_Updater(None, False, new_connector, fragment_id)
        updater.insert_bases(None, sequence)
        return updater.commit()


class Shared_Connector(object):
    """
    Represents a connection shared with another Connector object. Committing
    from a Shared_Connector object does not do anything. You must commit from
    the original Connector object.
    """

    def __init__(self, connector):
        self.__c = connector.conn
        self.__connector = connector

    @property
    def conn(self):
        """
        Returns sqlalchemy connection to database.
        """
        return self.__c

    # don't do anything on commit, we are a shared connector, so wait for parent
    # connector to commit.
    def commit(self):
        """
        No op
        """
        return self

    def new_connector(self, with_transaction):
        """
        Returns another Shared_Connector object sharing the same connection as
        the current Shared_Connector object.
        """
        # just reference arg to satisfy pylint
        if with_transaction is not None:
            with_transaction = None
        return Shared_Connector(self.__connector)
