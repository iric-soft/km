import os
import re
from itertools import groupby
from .Jellyfish import Jellyfish


def args_2_list_files(args):
    """Gather file names for ref. sequences."""
    if len(args) > 1:
        lst_files = args
    else:
        if os.path.isdir(args[0]):
            lst_files = [os.path.join(args[0], f) for f in os.listdir(args[0])]
        else:
            lst_files = args

    return(lst_files)


def target_2_seqfiles(target_fn):
    """Gather file names for ref. sequences."""
    return(args_2_list_files(target_fn))


def fasta_parser(fa_f):
    with open(fa_f, 'r') as f:
        groups = groupby(f, key=lambda x: x.startswith(">"))
        for is_header, lines in groups:
            if is_header:
                header = list(lines)[0].strip()
                sequence = "".join([l.strip() for l in next(groups)[1]])
                yield header, sequence


def file_2_seq(seq_f):
    sequences = []
    attributes = []
    for header, sequence in fasta_parser(seq_f):
        attr = {}
        for x in header.replace(">", "location=", 1).split("|"):
            k, v = x.split("=")
            attr[k.strip()] = v.strip()
        sequences.append(sequence.upper())
        attributes.append(attr)
    return sequences, attributes


def get_ref_kmer(ref_seq, ref_name, k_len):
    """ Load reference kmers. """

    ref_mer = []
    ref_set = set()
    for i in range(len(ref_seq) - k_len + 1):
        kmer = ref_seq[i:(i + k_len)]
        if kmer in ref_set:
            raise ValueError(
                "%s found multiple times in reference %s, at pos. %d" % (
                    kmer, ref_name, i)
            )
        ref_mer.append(kmer)
        ref_set.add(kmer)

    return ref_mer


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


def natsortkey(*args, rev_ix=[]):
    """Natural sorting of a string. For example: exon12 would
    come before exon2 with a regular sort, with natural sort
    the order would be exon2, exon12.
    """

    class reversor:
        def __init__(self, obj):
            self.obj = obj

        def __eq__(self, other):
            return other.obj == self.obj

        def __lt__(self, other):
            return self.obj > other.obj

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_split = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    split = lambda *l: tuple(alphanum_split(x) for x in l)
    reverse = lambda l, ix: tuple(reversor(x) if i in ix else x for i, x in enumerate(l))

    return reverse(split(*args), rev_ix)
