import re
import os
import os.path
import tempfile
import subprocess
from edge.models import Fragment, Genome
from edge.blast import BLAST_DB, genome_db_name
from edge.blast import Blast_Accession
from django.conf import settings
from django.core.management.base import BaseCommand


def fragment_fasta_fn(fragment):
    return '%s/edge-fragment-%s-nucl.fa' % (settings.NCBI_DATA_DIR, fragment.id)


def build_fragment_fasta(fragment):
    fn = fragment_fasta_fn(fragment)
    if not os.path.isfile(fn):  # have not built this fasta
        print 'building %s' % fn
        # this may take awhile, so do this first, so user interrupt does
        # not create an empty file
        sequence = fragment.indexed_fragment().sequence
        # be really lenient, convert any unknown bp to N
        sequence = re.sub(r'[^agctnAGCTN]', 'n', sequence)
        f = open(fn, 'w')
        f.write(">gnl|edge|%s %s\n%s\n" %
                (Blast_Accession.make(fragment), fragment.name, sequence))
        f.close()
    return fn


def build_db(fragments, dbname, refresh=True):
    if len(fragments) == 0:
        return

    if refresh is False and \
       (os.path.isfile(dbname+'.nal') or os.path.isfile(dbname+'.nsq')):
        print 'already built %s' % dbname
        return

    fns = []
    for fragment in fragments:
        fn = build_fragment_fasta(fragment)
        fns.append(fn)

    print 'concat fasta files for %s' % dbname
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        fafile = f.name
        for fn in fns:
            with open(fn) as inf:
                for line in inf:
                    f.write(line)

    print 'building blast db %s' % dbname
    cmd = "%s/makeblastdb -in %s -out %s " % (settings.NCBI_BIN_DIR, fafile, dbname)
    cmd += "-title edge -dbtype nucl -parse_seqids -input_type fasta"

    r = subprocess.check_output(cmd.split(' '))
    if 'Adding sequences from FASTA' not in r:
        print r

    os.unlink(fafile)


def build_genome_db(genome, refresh=False):
    fragments = list(genome.fragments.all())
    fragments.extend([op.fragment for op in genome.operations.all()])
    build_db(fragments, genome_db_name(genome), refresh=refresh)


def build_all_genome_dbs(refresh=False):
    for genome in Genome.objects.all():
        build_genome_db(genome, refresh=refresh)


def build_all_db():
    build_db(Fragment.objects.all(), BLAST_DB)
