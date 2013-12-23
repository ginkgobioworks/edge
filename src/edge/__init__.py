"""
Version control for genomic sequence changes. Allow derived genomes to inherit
annotations and changes from parent genomes.
"""

from edge.connector import Connector
from edge.genome import Genome
from edge.models import Annotation
from edge.fragment import Fragment, FragmentNotFound
from edge.io import IO
