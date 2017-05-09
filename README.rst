Edge
----

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
~~~~~~~~~~~~~~~~~~~

::

    make start-sim

Then check it out in your browser: http://localhost:9000/edge/#/genomes

Example: to import an genome, try

::

    make add-s288c-sim

The ``Makefile`` holds all the commands necessary for managing the server and
database, both in usage and development. Run ``make`` without arguments to see a list of commmands.

The Docker environment is defined in ``docker-compose.yml``. In that file, The ``edge-simple``
service builds a Docker image without NCBI and primer3. This takes less time and may be useful for
limited testing. Any make target invoked from the host machine that ends in ``-sim`` will be run
inside this environment.

To run all tests, and to use Edge for real, use the ``edge`` service. It works just like
``edge-simple``, except its commands end in ``-ext`` when run externally, e.g.

::

    make start-ext


Of course, any of these ``make`` targets can be run directly from a shell inside a container from
either service:

::

    you@localhost:edge$ docker-compose run --rm edge bash
    # Now you're inside the Docker container
    root@docker-image:/usr/src/edge# make test


Try it without Docker
~~~~~~~~~~~~~~~~~~~~~

Construct your virtual env and pip install dependencies (use
requirements.txt).

To start a server, first update ``src/server/settings.py`` to use either sqlite or MySQL. For MySQL,
create the appropriate databse. Then,

::

    python src/manage.py syncdb --noinput
    python src/manage.py migrate
    (cd example; gunzip ecoli-mg1655.gff.gz; gunzip yeast.gff.gz)
    python src/manage.py import_gff 'E. coli MG1655' example/ecoli-mg1655.gff
    python src/manage.py import_gff 'Saccharomyces cerevisiae' example/yeast.gff
    python src/manage.py runserver 0.0.0.0:8000

Then set your browser to http://localhost:8000/edge/

You can edit genome and fragment metadata, such as name, notes, circular
attributes, from the Django admin. Create a Django admin superuser, then set
your browser to http://localhost:8000/admin/

If you need BLAST support, run

::

    cd ncbi; ./install
    cd ../src; python manage.py build_edge_blastdb

If you need Primer3 support, run

::

    cd primer3; ./install


Development, testing, and release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When developing locally, you can run tests in the controlled environment of the docker container
from your local machine with ``make test-ext.``

Edge is versioned semantically. Continuous integration is done automatically on all branches through
Travis CI, and tagged commits are automatically released to PyPI. To release a new version,
bump the version number with the appropriate severity of the changes (major, minor, or patch), and
push the resulting tagged commits to the GitHub remote repo:

::
    you@localhost:edge$ make bump/patch-ext # Or bump/major, or bmp/minor
    you@localhost:edge$ git push --tags origin master

If you cannot push to master directly, do the same thing on a new branch and submit a pull request.
