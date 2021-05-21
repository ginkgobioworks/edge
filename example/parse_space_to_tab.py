# flake8: noqa
# Converting GFF format with space in lines starting with gi to tab

import sys
import re

fn = sys.argv[1]
f = open(fn, "r")

for l in f.read().split("\n"):
    if l.startswith("gi"):
        print re.sub(" ", "\t", l)
    else:
        print l
