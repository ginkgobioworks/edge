import os

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


def get_fragment_reference_fasta_gz_fn(fragment_id):
    return f"{settings.SEQUENCE_FILE_DIR}/edge-fragment-{fragment_id}.fa.gz"
