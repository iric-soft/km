import os
import sys
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


def file_2_seq(seq_f, k_len, individual=False):
    ref_seq = []
    ref_list = []
    ref_attr_list = []
    for line in open(seq_f, "r"):
        line = line.strip()
        if len(line) and line[0] == '>':
            if individual:
                if ref_seq:
                    ref_seq = "".join(ref_seq)
                    if len(ref_seq) >= k_len:
                        ref_list.append("".join(ref_seq))
                        ref_seq = []
                    else:
                        print >> sys.stderr, "Discarding", ref_attr_list.pop(-1)['exon']
                        ref_seq = []
                line = line.replace(">", "location=", 1)
                attr = {x.split("=")[0].strip():x.split("=")[1].strip() for x in line.split("|")}
                ref_attr_list.append(attr)
            continue
        ref_seq.append(line.upper())
    if individual:
        ref_list.append("".join(ref_seq))
        return ref_list, ref_attr_list
    else:
        ref_seq = ''.join(ref_seq)
        if len(ref_seq) >= k_len:
            return(ref_seq)
        else:
            raise IndexError(
                "Sequence in fasta is less than 31 nt"
            )


def get_ref_kmers(ref_seq, k_len, ref_name, check=True):
    """ Load reference kmers. """
    ref_mers = []
    ref_set = set()
    for i in range(len(ref_seq) - k_len + 1):
        kmer = ref_seq[i:(i + k_len)]
        if check and kmer in ref_set:
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


# Moved here to alleviate MutationFinder.py
def locate_reference(path, exon_seq, exon_names, exon_seq_ind, debug=False):
    exit_now = False
    import logging as log
    import numpy as np
    n_exons = len(exon_names)
    exon_num = range(n_exons)
    
    exon_table = np.empty([n_exons, len(path)])
    exon_table[:] = np.nan
    for i, k in enumerate(path):
        for e, seq_dict in exon_seq_ind.items():
            try:
                exon_table[e,i] = seq_dict[k]
            except KeyError:
                continue
    
    candidates = []
    for exon, seq, name, num in zip(exon_table, exon_seq, exon_names, exon_num):
        if 0 in exon:
            # Checks #1,2
            # Check 1: No more than 1 start;
            # Check 2: Array is sorted
            # Checks #1,2 fail in the case of ITDs; check 1 cannot fail without check 2 failing
            assert sum([x == 0 for x in exon]) < 2 or not np.all(np.diff(exon[~np.isnan(exon)]) > 0)
            
            start = np.nanargmin(exon)  # minimum is 0
            end = np.nanargmax(exon) # maximum is either length of exon or alternative splice site
            min_pos = np.nanmin(exon)
            max_pos = np.nanmax(exon)
            status = True if len(seq) == max_pos + 1 else False  # we got the whole sequence
            status_full = status and np.all(np.diff(exon[~np.isnan(exon)]) == 1)  # is not mutated 
            status_full = True if end + 1 == len(path) else status_full
            
            candidates.append(((num,), (name,), (start,), (end,), (status,), status_full))
            
            for exon2, seq2, name2, num2 in zip(exon_table, exon_seq, exon_names, exon_num):
                if np.any(np.isfinite(exon2)) and np.all(np.isnan(exon2[:end+21])):
                    # ^ if status is False, we don't know for how many steps the next exon was
                    #   wrongfully continued, so we set the maximum to 10 (31 - 21) (and that's generous)
                    # ^ same for True, because the next exon could be False
                    # Check #2 (again)
                    try:
                        assert np.all(np.diff(exon[~np.isnan(exon)]) > 0)
                    except AssertionError:
                        sys.stderr.write("Looks like we found an ITD in first exon in pair\n")
                    start2 = np.nanargmin(exon2)
                    end2 = np.nanargmax(exon2)
                    min_pos2 = np.nanmin(exon2)
                    max_pos2 = np.nanmax(exon2)
                    status2 = True if len(seq2) == max_pos2 - min_pos2 + 1 else False
                    status_full = True if status and status2 and end + 31 == start2 else False
                    candidates.append(((num, num2), (name, name2), (start, start2),
                                       (end, end2), (status, status2), status_full))
                    # Check #3
                    for exon3 in exon_table:
                        try:
                            assert not (np.any(np.isfinite(exon3)) and \
                                        np.all(np.isnan(exon3[:end2+31])))
                        except AssertionError:
                            sys.stderr.write("Found a 3rd exon, most probably nested, ignoring...\n")
    
    single_complete = [e for e in candidates if len(e[-2]) == 1 and e[-1]]
    double_all = [e for e in candidates if len(e[-2]) == 2]
    double_complete = []
    if double_all:
        if len(double_all) == 1:
            double_complete = double_all
        else:
            double_complete = [e for e in double_all if e[-1]]
            if not double_complete:
                double_complete = [e for e in double_all if sum(e[-2]) == 2 or sum(e[-2]) == 1]
                #^ decide later which references to keep (depending on the size of the ins/del)
                if not double_complete:
                    double_complete = double_all
    final_exons = single_complete + double_complete
    final_exon_ids = [e[0] for e in final_exons]
    final_exon_names = [e[1] for e in final_exons]
    
    if False:
        np.set_printoptions(threshold=np.nan)
        sys.stderr.write(str(exon_table[~np.all(np.isnan(exon_table), axis=1)]) + "\n")
        sys.stderr.write(
                "Exons: %s\n" % str(np.array(exon_names)[~np.all(np.isnan(exon_table), axis=1)])
                )
        sys.stderr.write("Length of path: " + str(len(path)) + "\n")
        for c in candidates:
            sys.stderr.write(str(c) + "\n")
        exit()
    
    return [e if len(e) == 2 else (e[0], None) for e in final_exon_ids]
