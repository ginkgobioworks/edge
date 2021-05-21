# flake8: noqa
# Converting GFF format that includes DNA as
#
# ##DNA chrI
# ##AGCTAGAT
# ##end-DNA
#
# to GFF with FASTA

import sys
import urllib
import re

fasta = []
curseq = []
started = False

fn = sys.argv[1]
f = open(fn, "r")


def unescape(s):
    return re.sub("\s", "_", urllib.unquote(s))


for l in f.read().split("\n"):
    l = l.strip()
    if l == "":
        continue
    if l.startswith("##"):
        if l.startswith("##DNA"):
            curseq.append(">%s" % unescape(l[6:]))
            started = True
        elif l.startswith("##end-DNA"):
            fasta.append("\n".join(curseq))
            curseq = []
            started = False
        else:
            if started:
                curseq.append(l[2:])
            else:
                print l
    else:
        if l.count("SGD\tchromosome") > 0:
            continue
        t = l.split("\t")
        t[0] = unescape(t[0])
        attrs = t[-1].split(";")
        for i, attr in enumerate(attrs):
            if len(attr.split("=")) != 2:
                attrs[i] = 'Unknown_Note="%s"' % attr.strip()
        t[-1] = ";".join(attrs)
        l = "\t".join(t)
        print l

print "##FASTA"
for f in fasta:
    print f
