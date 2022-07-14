import os


def make_required_dirs(path):
    dirn = os.path.dirname(path)
    try:
        original_umask = os.umask(0)
        os.makedirs(dirn, 0o777)
    except OSError:
        pass
    finally:
        os.umask(original_umask)
