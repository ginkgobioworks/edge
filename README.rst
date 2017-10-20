Edge
====

.. image:: https://travis-ci.org/ginkgobioworks/edge.svg?branch=master
    :target: https://travis-ci.org/ginkgobioworks/edge

Edge keeps structural changes between a genome and child genomes derived
from it. A user creates a modified genome by applying a sequence-based
operation, such as homologous recombination, to a parent genome. Users
can annotate or make corrections to sequences on a genome; Edge
automatically applies the changes to the appropriate regions on the
derived genomes. Edge does this efficiently: making a change on a parent
genome takes O(1) and is automatically propagated to the modified
genomes.

Edge uses O(D) amount of storage for each modified genome, where D is
the number of differences between a modified genome and its parent. The
current implementation additionally keeps a cache of annotations to base
pair numbers, but this cache is soft-data and is invalidated and
re-built on demand.

A modified genome can be re-created by re-applying operations to a new
genome (think git rebase). Currently, however, annotating a genome is
not an operation. Also, applying the same operation to a genome twice
results in a single child genome, not two.

Edge provides UIs to look at operations and changes and APIs for making
changes. Edge can export genome sequences and annotations as GFF files.
While Edge comes with a simple UI for browsing features and sequences,
the UI is primitive compared to other specialized applications.


Try it using Docker
-------------------
* Use ``docker-compose``:

To start edge server:
::
    
    docker-compose up
Then check it out in your browser: http://localhost:9000/edge/#/genomes
To import an genome using:
::
    docker-compose run edge python src/manage.py import_gff 'Saccharomyces cerevisiae' example/sc_s288c.gff
To run a shell inside the edge container
::
    docker-compose run --rm edge bash

* Alternatively, you can use the ``Makefile``:

To start edge server:
::

    make start-ext


To import an genome as an example:
::

    make add-s288c-ext


Or, if your app is already running, ``make add-s288c-ext_fast``, which will run the job without
trying to rebuild the image.

The ``Makefile`` holds all the commands necessary for managing the server and
database, both in usage and development. Run ``make`` without arguments to see a list of commmands.

The Docker environment is defined in ``docker-compose.yml``. Ese the ``edge`` service for your
commands.

Of course, any of these ``make`` targets can be run directly from a shell inside a container from
either service:

::

    you@localhost:edge$ docker-compose run --rm edge bash
    # Now you're inside the Docker container
    root@docker-image:/usr/src/edge# make test


Try it without Docker
---------------------

On your own machine, Construct your virtual environment and pip-install dependencies (use
``requirements.txt``).

To start a server, first update ``src/server/settings.py`` to use either sqlite or MySQL. For MySQL,
create the appropriate databse. Then,

::

    make migrate
    (cd example; gunzip ecoli-mg1655.gff.gz; gunzip yeast.gff.gz)
    python src/manage.py import_gff 'E. coli MG1655' example/ecoli-mg1655.gff
    python src/manage.py import_gff 'Saccharomyces cerevisiae' example/yeast.gff
    make run

Then set your browser to http://localhost:8000/edge/. Note the port is different than the Docker
case

If you need `NCBI BLAST`_ or Primer3_ support, you'll need to make sure the packages are installed
on your system. Debian and Ubuntu distributions provide binary versions of both of these packages.

Depending on where the NCBI BLAST tools and Primer3 are installed, you will probably need to tell
edge where to find them, using the following environment variables:

::

    NCBI_BIN_DIR       # Path to directory holding ncbi binaries, e.g. /usr/bin
    PRIMER3_BIN        # Path to primer3 binary, e.g. /usr/bin/primer3_core
    PRIMER3_CONFIG_DIR # Path to primer3 config directory, e.g. etc/primer3_config/


Then, to set up the edge BLAST db, from the ``src`` subdirectory,

::

    python manage.py build_edge_blastdb

.. _NCBI BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _Primer3: https://sourceforge.net/projects/primer3/


Editing data
------------

You can edit genome and fragment metadata, such as name, notes, circular
attributes, from the Django admin. Create a Django admin superuser, (see the ``superuser`` make
target), then set your browser to the ``/admin/`` endpoint of wherever you are running your dev
server.


Deploying to production
-----------------------

Do *not* use the Dockerfile as-is for production, or the ``make run`` task. Django's ``runserver``
command is not meant to run a production server. Instead, you'll need to spin up a production WSGI
server and run the Django projct with that, with your own settings. In this situation, it's better
to simply install the ``edge-genome`` python package on your deployed system and add it to your
deployment Django server's ``installed_apps`` setting.


Development, testing, and package release
-----------------------------------------

When developing locally, you can run tests in the controlled environment of the docker container
from your local machine with ``make test-ext.``

Edge is versioned semantically. Continuous integration is done automatically on all branches through
Travis CI, and tagged commits to master are automatically released to PyPI. To release a new version,
bump the version number with the appropriate severity of the changes (major, minor, or patch), and
push the resulting tagged commits to the GitHub remote repo:

::

    you@localhost:edge$ make bump/patch-ext # Or bump/major, or bmp/minor
    you@localhost:edge$ git push --tags origin master

If you cannot push to master directly, do the same thing on a new branch and submit a pull request.
