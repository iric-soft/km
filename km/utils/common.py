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

def file_2_fus_names(name_list):
    start_name = os.path.basename(name_list[0])
    end_name = os.path.basename(name_list[1])
    start_name = start_name[:-3]
    end_name = end_name[:-3]
    fus_names = []
    for i in range(len(file_2_exon_list(name_list[0]))):
        for j in range(len(file_2_exon_list(name_list[1]))):
            fus_names.append("%s_exon%d-%s_exon%d" % (start_name, i, end_name, j))
    return fus_names

def file_2_exon_list(seq_f):
    target_seq = []
    target_list = []
    first = True
    for line in open(seq_f, "r"):
        line = line.strip()
        if len(line) > 0 and line[0] == '>' and first:
            first = False
            continue;
        elif len(line) > 0 and line[0] == '>':
            target_seq = ''.join(target_seq)
            target_list.append(target_seq)
            target_seq = []
        else:
            target_seq.append(line)
    target_seq = ''.join(target_seq)
    target_list.append(target_seq)

    return target_list

def exons_2_fusion_seq(seq_list, reverse_compliment = True):
    start_list = file_2_exon_list(seq_list[0])
    end_list = file_2_exon_list(seq_list[1])
    fusion_seq = []
    for start in start_list:
        for end in end_list:
            if reverse_compliment:
                fusion_target = "%s" %start + "%s" % get_compliment(end)
            else:
                fusion_target = "%s" %start + "%s" % end
            fusion_seq.append(fusion_target)

    return fusion_seq

def get_compliment(seq):
    compliment = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    new_seq = []
    for base in seq:
        new_seq.append(compliment.get(base))
    new_seq = ''.join(new_seq)
    new_seq = new_seq[::-1]
    return new_seq
