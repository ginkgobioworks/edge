## Edge

Edge keeps structural changes between a genome and child genomes derived from it. A user creates a modified genome by applying a sequence based operation, such as homologous recombination, to a parent genome. Users can annotate or make corrections to sequences on a genome; Edge automatically applies the changes to the appropriate regions on the derived genomes. Edge does this efficiently: making a change on a parent genome takes O(1), and is automatically propagated to the modified genomes.

Edge uses O(D) amount of storage for each modified genome, where D is the number of differences between a modified genome and its parent. The current implementation additionally keeps a cache of annotations to base pair numbers, but this cache is soft-data and is invalidated and re-built on demand.

A modified genome can be re-created by re-applying operations to a new genome (think git rebase). Currently, however, annotating a genome is not an operation. Also, applying the same operation to a genome twice results in a single child genome, not two.

Edge provides UIs to look at operations and changes, and APIs for making changes. Edge can export genome sequences and annotations as GFF files.  While Edge comes with a simple UI for browsing features and sequences, the UI is primitive compared to other specialized applications. 

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

If you need Primer3 support, run

```
cd primer3; ./install
```
