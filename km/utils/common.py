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


def file_2_seq(seq_f, individual=False):
    ref_seq = []
    ref_list = []
    for line in open(seq_f, "r"):
        line = line.strip()
        if len(line) > 0 and line[0] == '>':
            if individual and ref_seq:
                ref_list.append("".join(ref_seq))
                ref_seq = []
                continue
            else:
                continue
        ref_seq.append(line)
    if individual:
        ref_list.append("".join(ref_seq))
        return(ref_list)
    else:
        ref_seq = ''.join(ref_seq)
        return(ref_seq)



def get_ref_kmers(ref_seq, k_len, ref_name):
    """ Load reference kmers. """
    ref_mers = []
    ref_set = set()
    for i in range(len(ref_seq) - k_len + 1):
        kmer = ref_seq[i:(i + k_len)]
        if kmer in ref_set:
            raise ValueError(
                "%s found multiple times in reference %s, at pos. %d" % (
                    kmer, ref_name, i)
            )
        
        ref_mers.append(kmer)
        ref_set.add(kmer)
    
    return(ref_mers)


def mean(v):
    if len(v) == 0:
        return 0
    else:
        return float(sum(v))/len(v)


def get_cov(db, ref_seq):
    jf = Jellyfish(db)
    cnt_stack = []

    i = 0
    count = 0
    cpt_count_0 = 0
    while i <= (len(ref_seq) - jf.k):
        kmer = ref_seq[i:i+jf.k]
        cnt = jf.query(kmer)
        cnt_stack += [int(cnt)]
        count += int(cnt)
        if int(cnt) == 0:
            cpt_count_0 += 1

        # print kmer, cnt
        i += 1

    return(count, len(ref_seq), min(cnt_stack), max(cnt_stack),
           mean(cnt_stack), len(cnt_stack), cpt_count_0)
