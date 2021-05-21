# flake8: noqa
import django

django.setup()

import sys
from edge import import_gff

import_gff(sys.argv[1], sys.argv[2])
