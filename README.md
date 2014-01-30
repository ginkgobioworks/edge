
Edge
----

Edge manages genomic changes and annotations. Edge allows a child genome to
derive from a parent genome by making a change to the parent genome. Children
genomes inherit annotations and changes from their parent genomes
appropriately.


Try it
------

Construct your virtual env, then try

```
python src/manage.py syncdb --noinput
python src/manage.py migrate
(cd example; gunzip ecoli-mg1656.gff.gz; gunzip yeast.gff.gz)
python src/manage.py import_gff 'E. coli MG1655' example/ecoli-mg1656.gff
python src/manage.py import_gff 'Saccharomyces cerevisiae' example/yeast.gff
python src/manage.py runserver 0.0.0.0:8000
```

Then set your browser to ```http://<your IP>:8000/edge/```

