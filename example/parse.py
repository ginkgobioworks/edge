# flake8: noqa
import sys
from BCBio import GFF


def parse(fn):
    in_handle = open(fn)
    for rec in GFF.parse(in_handle):
        for feature in rec.features:
            print feature
    in_handle.close()


if len(sys.argv) > 1:
    parse(sys.argv[1])
