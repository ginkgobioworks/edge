
## Edge

Edge keeps structural changes between a parent genome and its derived genomes.
Users create derived genomes by making sequence changes to the parent genome.
Derived genomes inherit annotations and fixes from their ancestors
automatically.


### Why?

Edge allows derived genomes to inherit annotations and changes to their
ancestors, even annotations and changes applied to an ancestor after the
creation of the derived genome. Edge does this efficiently: a change
(annotation or sequence change) to a parent genome takes O(1), and is
automatically propagated to genomes derived from the modified genome.

Edge uses O(D) amount of storage for each derived genome, where D is the number
of differences between the derived genome and its parent. The current
implementation additionally keeps a cache of annotations to base pair numbers,
but this cache is soft-data and is invalidated and re-built on demand.


### What Edge is NOT

Edge is not a genome sequence analyzer or viewer. There are much better tools
for those purposes. Edge's main goal is to keep track of genomic changes across
lineages in a structured fashion, so users can annotation, fix, and compare
genome sequences and features. While Edge comes with a simple viewer UI to look
at features and sequences, the UI is primitive. Edge does provide tools to
export genome sequences and annotations as GFF files.


### Try it

Construct your virtual env and pip install dependencies (use
requirements/{dev,core}.txt).

To start a server, first update src/server/settings.py to use either sqlite or
MySQL. For MySQL, create the appropriate databse. Then,

```
python src/manage.py syncdb --noinput
python src/manage.py migrate
(cd example; gunzip ecoli-mg1655.gff.gz; gunzip yeast.gff.gz)
python src/manage.py import_gff 'E. coli MG1655' example/ecoli-mg1655.gff
python src/manage.py import_gff 'Saccharomyces cerevisiae' example/yeast.gff
python src/manage.py runserver 0.0.0.0:8000
```

Then set your browser to ```http://<your IP>:8000/edge/```

You can edit genome and fragment metadata, such as name, notes, circular
attributes, from Django admin. Create a Django admin superuser, then set your
browser to ```http://<your IP>:8000/admin/```

If you need BLAST support, run

```
cd ncbi; ./install
cd ../src; python manage.py build_edge_blastdb
```
