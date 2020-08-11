# from a GFF and a FASTA file, create smaller GFFs, one for each sequence ID

import os
import sys

base_dir = sys.argv[1]
gff_fn = sys.argv[2]
fasta_fn = sys.argv[3]
gffs = []

# parse GFF and write out separate GFFs

cur_seqid = None
cur_seqid_fh = None

fh = open(gff_fn, "r")
for line in fh.readlines():
  if line.startswith("#"):
    continue
  tokens = line.split("\t")
  if len(tokens) < 5:
    continue
  seqid = tokens[0]

  if cur_seqid is None or seqid != cur_seqid: 
    if cur_seqid:
      cur_seqid_fh.close()
    print("gff %s" % cur_seqid)
    cur_seqid = seqid
    cur_seqid_fh = open("%s/%s.gff" % (base_dir, cur_seqid), "w")
    gffs.append(cur_seqid)

  cur_seqid_fh.write(line)

if cur_seqid:
  cur_seqid_fh.close()

fh.close()


# parse FASTA and write out separate GFFs

cur_seqid = None
cur_seqid_fh = None

fh = open(fasta_fn, "r")
for line in fh.readlines():
  if line.startswith(">"):
    tokens = line.split(" ")
    if len(tokens) < 1:
      continue
    seqid = tokens[0][1:]
    seqid = seqid.strip()
    if seqid not in gffs:
      print("%s not in GFF file" % seqid)
      continue

    if cur_seqid is None or seqid != cur_seqid: 
      if cur_seqid:
        cur_seqid_fh.close()
      print("fasta %s" % cur_seqid)
      cur_seqid = seqid
      cur_seqid_fh = open("%s/%s.fn" % (base_dir, cur_seqid), "w")

  cur_seqid_fh.write(line)

if cur_seqid:
  cur_seqid_fh.close()

fh.close()


# create combined gffs for importing

for gff in gffs:
  gff_fn = "%s/%s.gff" % (base_dir, gff)
  fna_fn = "%s/%s.fn" % (base_dir, gff)

  gff_fh = open(gff_fn, "r")
  fna_fh = open(fna_fn, "r")

  fh = open("%s/%s_combined.gff" % (base_dir, gff), "w")
  for line in gff_fh.readlines():
    fh.write(line)
  fh.write("##FASTA\n")
  for line in fna_fh.readlines():
    fh.write(line)

  gff_fh.close()
  fna_fh.close()
  fh.close()

  os.unlink(gff_fn)
  os.unlink(fna_fn)
