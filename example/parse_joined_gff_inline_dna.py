# Converting GFF format that includes DNA as
#
# ##DNA chrI
# ##AGCTAGAT
# ##end-DNA
#
# to GFF with FASTA

import sys

fasta = []
curseq = []
started = False

fn = sys.argv[1]
f = open(fn, 'r')

for l in f.read().split('\n'):
  if l.strip() == '':
    continue
  if l.startswith('##'):
    if l.startswith('##DNA'):
      curseq.append('>%s' % l[6:])
      started = True
    elif l.startswith('##end-DNA'):
      fasta.append('\n'.join(curseq));
      curseq = []
    else:
      if started:
        curseq.append(l[2:])
      else:
        print l
  else:
    if l.count('SGD\tchromosome') > 0:
      continue
    t = l.split('\t')
    attrs = t[-1].split(';')
    changed = False
    for i, attr in enumerate(attrs):
      if len(attr.split('=')) != 2:
        attrs[i] = 'Unknown_Note="%s"' % attr.strip()
        changed = True
    if changed:
      t[-1] = ';'.join(attrs)
      l = '\t'.join(t)
    print l

print '##FASTA'
for f in fasta:
  print f
