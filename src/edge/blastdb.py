import re
import os
import os.path
import tempfile
import subprocess
from edge.models import Fragment, Genome
from edge.blast import BLAST_DB, default_genome_db_name
from edge.blast import Blast_Accession
from django.conf import settings


def make_required_dirs(path):
    dirn = os.path.dirname(path)
    try:
        original_umask = os.umask(0)
        os.makedirs(dirn, 0o777)
    except OSError:
        pass
    finally:
        os.umask(original_umask)


def fragment_fasta_fn(fragment):
    return '%s/fragment/%s/%s/edge-fragment-%s-nucl.fa' % (settings.NCBI_DATA_DIR,
                                                           fragment.id % 1024,
                                                           (fragment.id >> 10) % 1024,
                                                           fragment.id)


def build_fragment_fasta(fragment, refresh=False):
    fn = fragment_fasta_fn(fragment)
    make_required_dirs(fn)

    if not os.path.isfile(fn) or refresh:  # have not built this fasta or need refresh
        print('building %s' % fn)
        # this may take awhile, so do this first, so user interrupt does
        # not create an empty file
        sequence = fragment.indexed_fragment().sequence
        # be really lenient, convert any unknown bp to N
        sequence = re.sub(r'[^agctnAGCTN]', 'n', sequence)
        if fragment.circular is True:
            sequence = sequence + sequence
        f = open(fn, 'w')
        f.write(">gnl|edge|%s %s\n%s\n" %
                (Blast_Accession.make(fragment), fragment.name, sequence))
        f.close()
    return fn


def build_db(fragments, dbname, refresh=True):
    if len(fragments) == 0:
        return None

    if refresh is False and \
       (os.path.isfile(dbname + '.nal') or os.path.isfile(dbname + '.nsq')):
        print('already built %s' % dbname)
        return dbname

    fns = []
    for fragment in fragments:
        fn = build_fragment_fasta(fragment, refresh)
        fns.append(fn)

    print('concat fasta files for %s' % dbname)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        fafile = f.name
        for fn in fns:
            with open(fn) as inf:
                for line in inf:
                    f.write(line)

    print('building blast db %s' % dbname)
    make_required_dirs(dbname)
    cmd = "%s/makeblastdb -in %s -out %s " % (settings.NCBI_BIN_DIR, fafile, dbname)
    cmd += "-title edge -dbtype nucl -parse_seqids -input_type fasta"

    r = subprocess.check_output(cmd.split(' '))
    if b'Adding sequences from FASTA' not in r:
        print(r)

    os.unlink(fafile)
    return dbname


def build_genome_db(genome, refresh=False):
    fragments = list(genome.fragments.all())
    dbname = build_db(fragments, default_genome_db_name(genome), refresh=refresh)
    genome.blastdb = dbname
    genome.save()


def check_and_build_genome_db(genome, refresh=False):
    if not genome.blastdb or refresh:
        build_genome_db(genome, refresh)


def build_all_genome_dbs(refresh=False):
    for genome in Genome.objects.all():
        build_genome_db(genome, refresh=refresh)


def build_all_db():
    build_db(Fragment.objects.all(), BLAST_DB)
