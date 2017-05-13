import math
from Bio.Seq import Seq

trans_table = 1       # standard translation table
min_protein_len = 100


def detect_orfs(seq):
    orf_list = []

    seq = Seq(seq)
    seq_len = len(seq)
    aa_len = int(math.floor(seq_len / 3.0))

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0

            # go through the translation and find end codons that follow a
            # start codon.
            while aa_start < trans_len and aa_start < aa_len:
                aa_end = trans.find("*", aa_start)
                has_stop = 1
                if aa_end == -1:
                    # no more stop codon, just abort...
                    break

                # we start looking for a M at the earliest at aa_end-aa_len+1,
                # since we don't want an ORF that's actually bigger than the
                # original sequence
                if aa_start < aa_end - aa_len + 1:
                    aa_start = aa_end - aa_len + 1
                start_codon = trans.find('M', aa_start, aa_end)

                # is there a start codon? and is it before end of sequence
                # (remember we doubled up the sequence earlier to detect orfs
                # crossing boundaries)
                if start_codon == -1 or start_codon >= aa_len:
                    assert(aa_end != -1)
                    aa_start = aa_end + 1
                    continue

                if aa_end - start_codon >= min_protein_len:
                    # the following start and end need to start with
                    # 1, not 0.
                    if strand == 1:
                        start = frame + start_codon * 3 + 1
                        end = frame + aa_end * 3 + has_stop * 3
                        if end > seq_len:
                            end = end % seq_len
                    else:
                        start = seq_len - frame - aa_end * 3 - has_stop * 3 + 1
                        end = seq_len - frame - start_codon * 3
                        if start < 0:
                            start = seq_len + start

                    f = dict(name='ORF frame ' + str(frame + 1),
                             start=start, end=end, strand=strand)
                    orf_list.append(f)

                aa_start = aa_end + 1

    return orf_list
