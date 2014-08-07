import re
import os
import os.path
import tempfile
import subprocess
from edge.models import Fragment
from edge.blast import BLAST_DB
from edge.blast import Blast_Accession
from Bio.Alphabet import IUPAC
from django.conf import settings
from django.core.management.base import BaseCommand


def fragment_fasta_fn(fragment):
    return '%s/edge-fragment-%s-nucl.fa' % (settings.NCBI_DATA_DIR, fragment.id)


def build_fragment_db():

    if Fragment.objects.count() == 0:
        return

    fns = []
    new = 0
    for fragment in Fragment.objects.all():
        fn = fragment_fasta_fn(fragment)
        fns.append(fn)
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
            new += 1

    print 'concat fasta files'
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        fafile = f.name
        for fn in fns:
            with open(fn) as inf:
                for line in inf:
                    f.write(line)

    print 'building blastdb'
    cmd = "%s/makeblastdb -in %s -out %s " % (settings.NCBI_BIN_DIR, fafile, BLAST_DB)
    cmd += "-title edge -dbtype nucl -parse_seqids -input_type fasta"

    r = subprocess.check_output(cmd.split(' '))
    if 'Adding sequences from FASTA' not in r:
        print r

    os.unlink(fafile)


class Command(BaseCommand):
    def handle(self, *args, **options):
        build_fragment_db()
