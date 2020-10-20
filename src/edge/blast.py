import os
import tempfile
import subprocess
from functools import lru_cache
from django.conf import settings
from edge.models import Fragment
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML


BLAST_DB = "%s/edge-nucl" % settings.NCBI_DATA_DIR


def default_genome_db_name(genome):
    return "%s/genome/%s/%s/edge-genome-%d-nucl" % (settings.NCBI_DATA_DIR,
                                                    genome.id % 1024, (genome.id >> 10) % 1024,
                                                    genome.id)


class Blast_Accession(object):

    @staticmethod
    def make(fragment):
        return '%d/%d' % (fragment.id, fragment.indexed_fragment().length)

    def __init__(self, accession):
        v = str(accession).split('/')
        self.fragment_id = int(v[0])
        if len(v) > 1:
            self.fragment_length = int(v[1])
        else:
            self.fragment_length = None

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
        return self.identities() * 1.0 / self.alignment_length()


def inverse_match(m):
    return ''.join([' ' if x == '|' else 'X' for x in m])


@lru_cache(maxsize=200)
def blast(dbname, blast_program, query):

    infile = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        infile = f.name
        f.write(">Query\n%s\n" % query)

    outfile = "%s.out.xml" % infile
    if blast_program == 'tblastn':
        blast_cl = NcbitblastnCommandline(query=infile, db=dbname,
                                          word_size=6, outfmt=5, out=outfile)
    else:
        blast_cl = NcbiblastnCommandline(query=infile, db=dbname,
                                         word_size=6, outfmt=5, out=outfile)

    cl = str(blast_cl)
    cl = "%s/%s" % (settings.NCBI_BIN_DIR, cl)
    r = subprocess.call(cl.split(" "))
    os.unlink(infile)

    if r != 0:
        print("Blast failed: %s" % cl)
        return []

    results = []
    with open(outfile, "r") as f:
        blast_record = NCBIXML.read(f)
        for alignment in blast_record.alignments:
            accession = Blast_Accession(alignment.accession)
            for hsp in alignment.hsps:
                if accession.fragment_length is not None:
                    if hsp.sbjct_start > accession.fragment_length and \
                       hsp.sbjct_end > accession.fragment_length:
                        continue
                    # don't apply '% accession.fragment_length' to
                    # sbjct_start/end. Blast_Result#strand compares sbjct_start
                    # and sbjct_end to determine which strand the hit is on.
                    # Caller should just handle when sbjct_start/end is greater
                    # than fragment length. alternatively, we can store strand
                    # explicit, but that also creates complexity when using
                    # sbjct_start/end coordinates.

                f = Blast_Result(fragment_id=accession.fragment_id,
                                 fragment_length=accession.fragment_length,
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


def blast_genome(genome, blast_program, query):
    dbname = genome.blastdb
    if not dbname:
        return []
    results = blast(dbname, blast_program, query)

    genome_fragment_ids = [f.id for f in genome.fragments.all()]
    return [r for r in results if r.fragment_id in genome_fragment_ids]
