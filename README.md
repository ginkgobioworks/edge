
## Edge

Edge keeps track of proposed structural changes between a parent genome and its
derived genomes. Users can create a derived genome by making changes to the
sequence of the parent genome. Derived genomes inherit annotations from parents
whenever appropriate.


### Why?

Edge efficiently keeps track of structural changes between genomes, and at the
same time allows derived genomes to inherit annotations from parents, even
annotations applied to a parent after the derived genome was created.

The necessary storage cost to store a derived genome is O(D), where D is the
number of differences between the derived genome and its parent. The current
implementation additionally uses a cache of annotation to base pair number for
each genome. While this cache is O(N), where N is the number of annotations, it
is soft-data and can be discarded.


### Try it

Construct your virtual env and pip install dependencies (use
requirements/{dev,core}.txt).

To start a server,

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

