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


def genome_db_name(genome):
    return "%s/edge-genome-%d-nucl" % (settings.NCBI_DATA_DIR, genome.id)


class Blast_Accession(object):

    @staticmethod
    def make(fragment):
        return '%s' % fragment.id

    def __init__(self, accession):
        self.fragment_id = int(accession)

    @property
    def fragment(self):
        return Fragment.objects.get(pk=self.fragment_id)


class Blast_Result(object):

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def to_dict(self):
        return self.__dict__

    @property
    def fragment(self):
        return Fragment.objects.get(pk=self.fragment_id)

    def strand(self):
        start = self.subject_start
        end = self.subject_end
        if start < end:
            return 1
        else:
            return -1

    def alignment_length(self):
        return len(self.alignment['match'])

    def identities(self):
        identities = 0
        for q, s in zip(self.alignment['query'], self.alignment['subject']):
            if q.lower() == s.lower():
                identities += 1
        return identities

    def identity_ratio(self):
        return self.identities()*1.0/self.alignment_length()


def inverse_match(m):
    return ''.join([' ' if x == '|' else 'X' for x in m])


def blast(dbname, blast_program, query, evalue_threshold=0.001):

    infile = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        infile = f.name
        f.write(">Query\n%s\n" % query)

    outfile = "%s.out.xml" % infile
    if blast_program == 'tblastn':
        blast_cl = NcbitblastnCommandline(query=infile, db=dbname,
                                          evalue=evalue_threshold,
                                          word_size=6, outfmt=5, out=outfile)
    else:
        blast_cl = NcbiblastnCommandline(query=infile, db=dbname,
                                         evalue=evalue_threshold,
                                         word_size=6, outfmt=5, out=outfile)

    cl = str(blast_cl)
    cl = "%s/%s" % (settings.NCBI_BIN_DIR, cl)
    # print cl
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
                f = Blast_Result(fragment_id=accession.fragment_id,
                                 hit_def=alignment.hit_def,
                                 query_start=hsp.query_start,
                                 query_end=hsp.query_end,
                                 subject_start=hsp.sbjct_start,
                                 subject_end=hsp.sbjct_end,
                                 evalue=hsp.expect,
                                 alignment=dict(query=hsp.query,
                                                match=hsp.match,
                                                matchi=inverse_match(hsp.match),
                                                subject=hsp.sbjct))
                results.append(f)

    os.unlink(outfile)
    return results


def blast_genome(genome, blast_program, query, evalue_threshold=0.001):
    results = blast(genome_db_name(genome), blast_program, query,
                    evalue_threshold=evalue_threshold)
    genome_fragment_ids = [f.id for f in genome.fragments.all()]
    genome_fragment_ids.extend([op.fragment_id for op in genome.operations.all()])
    return [r for r in results if r.fragment_id in genome_fragment_ids]
