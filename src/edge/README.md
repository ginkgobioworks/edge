We keep edges between chunks. Each edge is labeled with a fragment ID. A
fragment F has two consecutive chunks A and B if there is an edge between A and
B labeled with F, or if there is an edge between A and B labeled with one of
F's ancestors and there is no edge between A and any other chunk labeled with a
more immediate ancestor of F.

Using these edges, you can walk through all the chunks of a fragment, w/o the
need to keep location indices. The number of edges you need for a lineage of
genomes is the number of annotations on base genome plus number of genomic
changes and new annotations on children of the base genome, whereas the number
of location indices you need would equal to total number of inherited and new
annotations on all genomes.
