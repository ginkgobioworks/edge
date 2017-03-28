import re
import sys

f = open(sys.argv[1],"r")
lines = f.readlines()
f.close()

print "##GFF"
for line in lines:
  words = line.strip().split("\t")
  type = words[1]
  name = words[3]
  alias = words[4]
  chr = int(words[8])
  start = int(words[9])
  end = int(words[10])
  notes = ""
  if len(words) >= 16:
    notes = words[15].rstrip("/")
    notes = re.sub("\s", "%20", notes)

  gene = ""
  if type == "ORF":
    type = "gene"
    gene = "gene=%s;" % alias

  strand = "+"
  if start > end:
    strand = "-"
    x = end
    end = start
    start = x

  print "chr%s\tAmyris\t%s\t%s\t%s\t.\t%s\t0\tID=%s;Name=%s;%sNote=%s" % (chr, type, start, end, strand, name, alias, gene, notes)

f = open(sys.argv[2],"r")
lines = f.readlines()
f.close()

print "##FASTA"
for line in lines:
  line = line.strip()
  if line[0] == ">":
    print ">chr%s" % line[1:]
  else:
    print line
