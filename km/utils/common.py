import os
from Jellyfish import Jellyfish


def args_2_list_files(args):
    """Gather file names for ref. sequences."""
    if len(args) > 1:
        lst_files = args
    else:
        if os.path.isdir(args[0]):
            lst_files = map(
                lambda f: os.path.join(args[0], f),
                os.listdir(args[0]))
        else:
            lst_files = args

    return(lst_files)


def target_2_seqfiles(target_fn):
    """Gather file names for ref. sequences."""
    return(args_2_list_files(target_fn))


def file_2_seq(seq_f):
    target_seq = []
    for line in open(seq_f, "r"):
        line = line.strip()
        if len(line) > 0 and line[0] == '>':
            continue
        target_seq.append(line)
    target_seq = ''.join(target_seq)

    return(target_seq)


def get_target_kmers(target_seq, k_len, target_name):
    """ Load reference kmers. """
    target_kmers = []
    for i in range(len(target_seq) - k_len + 1):
        kmer = target_seq[i:(i + k_len)]
        if kmer in target_kmers:
            raise ValueError(
                "%s found multiple times in reference %s, at pos. %d" % (
                    kmer, target_name, i)
            )

        target_kmers.append(kmer)

    return(target_kmers)


def mean(v):
    if len(v) == 0:
        return 0
    else:
        return float(sum(v))/len(v)


def get_cov(db, target_seq):
    jf = Jellyfish(db)
    cnt_stack = []

    i = 0
    count = 0
    cpt_count_0 = 0
    while i <= (len(target_seq) - jf.k):
        kmer = target_seq[i:i+jf.k]
        cnt = jf.query(kmer)
        cnt_stack += [int(cnt)]
        count += int(cnt)
        if int(cnt) == 0:
            cpt_count_0 += 1

        # print kmer, cnt
        i += 1

    return(count, len(target_seq), min(cnt_stack), max(cnt_stack),
            mean(cnt_stack), len(cnt_stack), cpt_count_0)
