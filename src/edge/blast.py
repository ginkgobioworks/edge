import os
import re
import tempfile
import subprocess
from django.conf import settings
from edge.models import Fragment
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML


BLAST_DB = "%s/edge-nucl" % settings.NCBI_DATA_DIR


class Blast_Accession(object):

    @staticmethod
    def make(fragment):
        return '%s' % fragment.id

    def __init__(self, accession):
        self.fragment_id = int(accession)

    @property
    def fragment(self):
        return Fragment.objects.get(pk=self.fragment_id)


def blast(query, blast_program, evalue_threshold=0.001):

    infile = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        infile = f.name 
        f.write(">Query\n%s\n" % query)

    outfile = "%s.out.xml" % infile
    if blast_program == 'tblastn':
        blast_cl = NcbitblastnCommandline(query=infile, db=BLAST_DB,
                                          evalue=evalue_threshold, word_size=6, outfmt=5, out=outfile)
    else:
        blast_cl = NcbiblastnCommandline(query=infile, db=BLAST_DB,
                                         evalue=evalue_threshold, word_size=6, outfmt=5, out=outfile)

    cl = str(blast_cl)
    cl = "%s/%s" % (settings.NCBI_BIN_DIR, cl)
    print cl
    r = subprocess.call(cl.split(" "))
    os.unlink(infile)

    if r != 0:
        print "Blast failed: %s" % cl
        return []

    results = []
    with open(outfile, "r") as f:
        blast_record = NCBIXML.read(f)
        for alignment in blast_record.alignments:
            accession = Blast_Accession(alignment.accession)
            for hsp in alignment.hsps:
                f = dict(fragment_id=accession.fragment_id,
                         hit_def=alignment.hit_def,
                         query_start=hsp.query_start,
                         query_end=hsp.query_end,
                         subject_start=hsp.sbjct_start,
                         subject_end=hsp.sbjct_end,
                         evalue=hsp.expect,
                         alignment=dict(query=hsp.query, match=hsp.match, subject=hsp.sbjct))
                results.append(f)

    os.unlink(outfile)
    return results
