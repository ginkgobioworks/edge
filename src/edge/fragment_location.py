import pickle
#import pandas as pd
from sortedcontainers import SortedDict

from edge.models.chunk import Chunk


def unpickle_fl_file(fn):
    """
    Unpickles FL file into SortedDict.

    Parameters
    ----------
    fn: string
        The pickled filename to be read

    Returns
    -------
    SortedDict
        {Start Coordinate: Chunk ID, ...}
    """
    with open(fn, 'rb') as f:
        coord_to_chunk_id = pickle.load(f)
    return coord_to_chunk_id

def pickle_fl_file(fn, coord_to_chunk_id):
    """
    Pickles SortedDict into FL file.

    Parameters
    ----------
    fn: string
        The pickled filename to be written to
    coord_to_chunk_id: SortedDict{int: int}
        {Start Coordinate: Chunk ID, ...}

    Returns
    -------
    None
    """
    with open(fn, 'wb') as f:
        pickle.dump(coord_to_chunk_id, f)

def new_fl_file(fragment_id, chunk_ids_and_sizes, dirn='.'):
    """
    Creates and pickles empty FL file.

    Parameters
    ----------
    fragment_id: int
        ID of the associated fragment for file naming
    chunk_ids_and_sizes: list[tuple(int, int)]
        [(1st Chunk ID, 1st Chunk Size), ...]
    dirn: string (optional)
        Directory to store file in

    Returns
    -------
    string
        The pickled FL filename
    """
    coord_to_chunk_id = SortedDict()
    coord = 1
    for chunk_id, size in chunk_ids_and_sizes:
        coord_to_chunk_id[coord] = chunk_id
        coord += size
    coord_to_chunk_id[coord] = None

    fn = "%s/%s.p" % (dirn, fragment_id)
    pickle_fl_file(fn, coord_to_chunk_id)
    return fn

def reconstruct_fragment_from_fl_file(fn):
    """
    Reconstructs fragment sequence from FL file.

    Parameters
    ----------
    fn: string
        The pickled FL filename

    Returns
    -------
    string
        The fragment's sequence
    """
    coord_to_chunk_id = unpickle_fl_file(fn)

    sequences_list = []
    for chunk_id in coord_to_chunk_id.values():
        if chunk_id is None:
            break
        chunk = Chunk.objects.get(pk=chunk_id)
        sequences_list.append(chunk.get_sequence())

    return ''.join(sequences_list)

def reconstruct_fragment_from_fl_file_by_coordinate(fn, start, end):
    """
    Reconstructs fragment sequence from FL file based on given coordinates.

    Parameters
    ----------
    fn: string
        The pickled FL filename
    start: int
        The 1-indexed start coordinate
    end: int
        The 1-indexed end coordinate

    Returns
    -------
    string
        The fragment's sequence
    """
    coord_to_chunk_id = unpickle_fl_file(fn)

    # force inputs to fit bounds
    start = max(1, start)
    end = min(coord_to_chunk_id.keys()[-1] - 1, end)

    # find indices and positions of chunks involved
    start_chunk_index = coord_to_chunk_id.bisect_right(start) - 1
    start_chunk_pos = coord_to_chunk_id.keys()[start_chunk_index]
    end_chunk_index = coord_to_chunk_id.bisect_right(end) - 1
    end_chunk_pos = coord_to_chunk_id.keys()[end_chunk_index]

    sequences_list = []
    for k in coord_to_chunk_id.keys()[start_chunk_index: end_chunk_index + 1]:
        # get chunk and sequence
        chunk_id = coord_to_chunk_id[k]
        if chunk_id is None:
            break
        chunk = Chunk.objects.get(pk=chunk_id)
        seq = chunk.get_sequence()

        # slice sequences (if necessary) for first and last chunks
        if k == start_chunk_pos:
            seq = seq[(start - k):]
        elif k == end_chunk_pos:
            seq = seq[0:(end - k + 1)]
        sequences_list.append(seq)

    return ''.join(sequences_list)

"""def run_baselines_on_gtf(fn):
    def estimate_size_by_contig(contig, fn):
        start_time = time.time()
        column_names = ['contig', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 'extra']
        df = pd.read_csv("../../hg38.knownGene.gtf", delimiter='\t', header=None, names=column_names)
        df = df[df['contig'] == contig]
        print("--------------------------")
        print(f"Number of {contig} associated rows:", len(df.index))

        pivots = set(pd.unique(df[['start', 'end']].values.ravel('K')))
        pivots.add(1)
        n_chunks = len(pivots)
        print("Number of chunks:", n_chunks, "\n")

        construction_start_time = time.time()
        coord_to_chunk_id = SortedDict()
        for start_coord in pivots:
            coord_to_chunk_id[start_coord] = start_coord + 1

        fn = "%s/%s" % (".", "test")
        pickle_start_time = time.time()
        pickle_fl_file(fn, coord_to_chunk_id)
        print("Time to pandas parse:", construction_start_time - start_time, "seconds")
        print("Time to construct dictionaries:", pickle_start_time - construction_start_time, "seconds")
        print("Time to pickle:", time.time() - pickle_start_time, "seconds")
        print("Time for all:", time.time() - start_time, "seconds")
        print("Time per chunk", (time.time() - start_time) / n_chunks, "seconds\n")
        
        file_size = os.path.getsize(fn)

        print("Size of file:", file_size / (10**6), "MB")
        print("Bytes per chunk:", file_size / n_chunks)

        return n_chunks, file_size

    total_start_time = time.time()
    total_chunks = 0
    total_bytes = 0
    for i in range(1, 23):
        contig = f"chr{i}"
        new_chunks, new_bytes = estimate_size_by_contig(contig, fn)
        total_chunks += new_chunks
        total_bytes += new_bytes

    print("--------------------------")
    print("Over 22 chromosomes:")
    print("Total chunks:", total_chunks)
    print("Total bytes:", total_bytes)
    print("Total bytes per chunk", total_bytes / total_chunks)
    print("Total time:", time.time() - total_start_time, "seconds")
    print("Total time per chunk", (time.time() - total_start_time) / total_chunks, "seconds")
"""