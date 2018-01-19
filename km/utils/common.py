import os


def target_2_seqfiles(target_fn):
    """Gather file names for ref. sequences."""

    if len(target_fn) > 1:
        seq_files = target_fn
    else:
        if os.path.isdir(target_fn[0]):
            seq_files = map(
                lambda f: os.path.join(target_fn[0], f),
                os.listdir(target_fn[0]))
        else:
            seq_files = target_fn

    return(seq_files)


def file_2_seq(seq_f):
    ref_seq = []
    for line in open(seq_f, "r"):
        line = line.strip()
        if line[0] == '>':
            continue
        ref_seq.append(line)
    ref_seq = ''.join(ref_seq)

    return(ref_seq)


def get_ref_kmer(ref_seq, k_len, ref_name):
    """ Load reference kmers. """
    ref_mer = []
    for i in range(len(ref_seq) - k_len + 1):
        kmer = ref_seq[i:(i + k_len)]
        if kmer in ref_mer:
            raise ValueError(
                "%s found multiple times in reference %s, at pos. %d" % (
                    kmer, ref_name, i)
            )

        ref_mer.append(kmer)

    return(ref_mer)
