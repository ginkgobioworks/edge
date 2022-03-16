import pickle
import os, time
import pandas as pd
from sortedcontainers import SortedDict

from edge.models.chunk import Chunk


def unpickle_fl_file(fn):
    """
    Unpickles FL file into objects.

    Parameters
    ----------
    fn: string
        The pickled filename to be read

    Returns
    -------
    SortedDict
        {Start Coordinate: Chunk ID, ...}
    dict
        {Chunk ID: Start Coordinate, ...}
    """
    with open(fn, 'rb') as f:
        coord_to_chunk_id = pickle.load(f)
        chunk_id_to_coord = pickle.load(f)
    return coord_to_chunk_id, chunk_id_to_coord

def pickle_fl_file(fn, coord_to_chunk_id, chunk_id_to_coord):
    """
    Pickles objects into FL file.

    Parameters
    ----------
    fn: string
        The pickled filename to be written to
    coord_to_chunk_id: SortedDict{int: int}
        {Start Coordinate: Chunk ID, ...}
    chunk_id_to_coord: dict{int: int}
        {Chunk ID: Start Coordinate, ...}

    Returns
    -------
    None
    """
    with open(fn, 'wb') as f:
        pickle.dump(coord_to_chunk_id, f)
        pickle.dump(chunk_id_to_coord, f)

def new_fl_file(fragment_id, chunk_ids_and_sizes, dir='.'):
    """
    Creates and pickles empty FL file.

    Parameters
    ----------
    fragment_id: int
        ID of the associated fragment for file naming
    chunk_ids_and_sizes: list[tuple(int, int)]
        [(1st Chunk ID, 1st Chunk Size), ...]
    dir: string (optional)
        Directory to store file in

    Returns
    -------
    string
        The pickled FL filename
    """
    coord_to_chunk_id = SortedDict()
    chunk_id_to_coord = {}
    coord = 1
    for chunk_id, size in chunk_ids_and_sizes:
        coord_to_chunk_id[coord] = chunk_id
        chunk_id_to_coord[chunk_id] = coord
        coord += size
    coord_to_chunk_id[coord] = None

    fn = "%s/%s" % (dir, fragment_id)
    pickle_fl_file(fn, coord_to_chunk_id, chunk_id_to_coord)
    return fn

def update_fl_file(fn, old_chunk_id_to_new_chunk_ids_and_sizes):
    """
    Updates existing FL file.

    Parameters
    ----------
    fn: string
        The pickled FL filename
    old_chunk_id_to_new_chunk_ids_and_sizes: dict{int: list[tuple(int, int)]}
        {Old Chunk ID: [(1st New Chunk ID, 1st New Chunk Size), ...]}

    Returns
    -------
    None
    """
    coord_to_chunk_id, chunk_id_to_coord = unpickle_fl_file(fn)

    for old_chunk_id, new_chunk_ids_and_sizes in old_chunk_id_to_new_chunk_ids_and_sizes:
        coord = chunk_id_to_coord.pop(old_chunk_id, None)
        if coord is None:
            continue
        stored_old_chunk_id = coord_to_chunk_id.pop(coord)
        assert old_chunk_id == stored_old_chunk_id

        for chunk_id, size in new_chunk_ids_and_sizes:
            coord_to_chunk_id[coord] = chunk_id
            chunk_id_to_coord[chunk_id] = coord
            coord += size

    pickle_fl_file(fn, coord_to_chunk_id, chunk_id_to_coord)

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
    coord_to_chunk_id, _ = unpickle_fl_file(fn)

    sequences_list = []
    for chunk_id in coord_to_chunk_id.values():
        chunk = Chunk.objects.get(pk=chunk_id)
        sequences_list.append(chunk.get_sequence())

    return ''.join(sequences_list)

def run_baselines_on_gtf(fn):
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
        chunk_id_to_coord = {}
        for start_coord in pivots:
            coord_to_chunk_id[start_coord] = start_coord + 1
            chunk_id_to_coord[start_coord + 1] = start_coord + 1

        fn = "%s/%s" % (".", "test")
        pickle_start_time = time.time()
        pickle_fl_file(fn, coord_to_chunk_id, chunk_id_to_coord)
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
