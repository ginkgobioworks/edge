# takes two sequences, strips whitespaces away, computes reverse
# complement sequence of second primer

a = "A T G A C A A G C G A A C C A G A G T T T C A G C A G G C T T A C G A T G A G A T C G T T T C T T C T G T G G A G G A T T C C A A A A T T T T T G A A A A A T T C C C A C A G T A T A A A A A A G T G T T A C C T A T T G T T T C T G T C C C G G A G A G G A T C A T T C A A T T C A G G G T C A C G T G G G A A A A T G A T A"

b = "T A C A A A C A C C T T G C C A T C A T T G G T C A A G G G G G C C A A C A T T G C C A G C T T C G T C A T G G T G G C T G A C G C A A T G C T T G A C C A G G G A G A C G T T T T T T A G C C G T A A G C G C T A T T T T C T T T T T G T T C G T A A C T A T C T G T G T A T G T A G T A G T G T A A T C T A C T T T T A A T"

import re
from Bio.Seq import Seq

a = re.sub(r'\s', '', a)
b = re.sub(r'\s', '', b)

print a
print str(Seq(b).reverse_complement())
