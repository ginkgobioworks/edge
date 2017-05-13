import sys

gff = []
fasta = []

for fn in sys.argv[1:]:
  f = open(fn + '.gff', 'r')
  g = []
  s = []
  fasta_mode = False
  for l in f.read().split('\n'):
    if l.startswith('##FASTA'):
      fasta_mode = True
    elif len(l) == 0 or l[0] == '#':
      continue
    elif 'annotation\tremark' in l:
      continue
    elif fasta_mode is False:
      t = l.split('\t')
      t[0] = fn
      l = '\t'.join(t)
      g.append(l)
    else:
      if l[0] == '>':
        l = '>%s %s' % (fn, l[1:])
      s.append(l)
  gff.append('\n'.join(g))
  fasta.append('\n'.join(s))

for g in gff:
  print g
print '##FASTA'
for f in fasta:
  print f
