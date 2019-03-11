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
    ref_attr_list = []
    for line in open(seq_f, "r"):
        line = line.strip()
        if len(line) and line[0] == '>':
            if individual:
                line = line.replace(">", "location=", 1)
                attr = {x.split("=")[0].strip():x.split("=")[1].strip() for x in line.split("|")}
                ref_attr_list.append(attr)
                if ref_seq:
                    ref_list.append("".join(ref_seq))
                    ref_seq = []
            continue
        ref_seq.append(line.upper())
    if individual:
        ref_list.append("".join(ref_seq))
        return ref_list, ref_attr_list
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


# Moved here to alleviate MutationFinder.py
def locate_reference(path, exon_seq, exon_names, exon_seq_ind, debug=False):
    exit_now = False
    import sys
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
            # Checks #1 and #2 fail in the case of ITDs, and check 1 cannot fail without check 2 failing
            # Check #1
            check_1 = True
            try:
                assert sum([x == 0 for x in exon]) < 2  # no more than 1 start
            except AssertionError:
                check_1 = False
            # Check #2
            check_2 = True
            try:
                assert np.all(np.diff(exon[~np.isnan(exon)]) > 0)  # array is sorted
            except AssertionError:
                #sys.stderr.write("ERROR: Check 2 fail\n" + str(exon) + "\n" + str(path) + "\n" +
                #        str(len(seq)) + "\n" + str(name) + "\n")
                check_2 = False
            # Check #1-2
            assert check_1 or (not check_1 and not check_2)
            
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
                    if not check_2:
                        sys.stderr.write("Looks like we found an ITD in first exon in pair\n")
                        #exit_now = True
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
                            #sys.stderr.write("ERROR: Check 3 fail\n" + str(exon3) + "\n" +
                            #        "\n" + str(len(seq)) + "\n" + str(name) + "\n")
                            sys.stderr.write("Found a 3rd exon, most probably nested, ignoring...\n")
                            #exit_now = True
    
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
                # decide later which references to keep (depending on the length of the ins/del)
                if not double_complete:
                    double_complete = double_all
    final_exons = single_complete + double_complete
    final_exon_ids = [e[0] for e in final_exons]
    final_exon_names = [e[1] for e in final_exons]
    
    if exit_now:
        np.set_printoptions(threshold=np.nan)
        #print exon_table
        #print exon_table[~np.all(np.isnan(exon_table), axis=1)]
        #print np.array(exon_names)[~np.all(np.isnan(exon_table), axis=1)]
        sys.stderr.write(str(len(path)) + "\n")
        for c in candidates:
            sys.stderr.write(str(c) + "\n")
        #exit()
    
    return [e if len(e) == 2 else (e[0], None) for e in final_exon_ids]
    #
    ##debug = True
    #start_kmers = exon_seq_ind
    #def find_all(path, exons_in_order, exon_alternatives):
    #    if not path:
    #        exon_alternatives.append(exons_in_order)
    #        if debug:
    #            sys.stderr.write("RETURNING {}\n".format(exon_alternatives))
    #        return exon_alternatives
    #    exons_to_explore = []
    #    for i, k in enumerate(path):
    #        if len(path) == 1:  # specific scenario
    #            i += 1
    #            assert k in start_kmers
    #            if debug:
    #                sys.stderr.write("START " + str(k) + "\n")
    #            l = len(exons_to_explore)
    #            for e in exon_num:
    #                if exon_seq[e][0] == k:
    #                    exons_to_explore.append([e, len(exon_seq[e]), 1, 0, False, l])
    #        if exons_to_explore:
    #            if debug:
    #                sys.stderr.write(str(exons_to_explore) + "\n")
    #            for j, ee in enumerate(exons_to_explore):
    #                cur_seq = exon_seq[ee[0]]
    #                if ee[3] < 0:
    #                    ee[3] -= 1
    #                elif ee[4] or ee[1] == ee[2]:
    #                    ee[2] += 1
    #                    ee[4] = True
    #                elif len(cur_seq) > ee[2] + ee[3] and k == cur_seq[ee[2] + ee[3]]:
    #                    ee[2] += 1
    #                    ee[2] += ee[3] # we had a snp or mnp
    #                    ee[3] = 0
    #                    if ee[1] == ee[2]:
    #                        ee[4] = True
    #                    else:
    #                        ee[4] = False
    #                elif k == exon_seq[ee[0]][-1]:  # we had an indel/ins/del
    #                    ee[2] += ee[3]
    #                    ee[3] = 0
    #                    ee[4] = True
    #                else:  # we might have an indel/ins/del
    #                    if ee[1] - ee[2] + 1 <= 31:  # No time to recover from a point mutation
    #                        ee[3] = -1
    #                        ee[4] = None
    #                    else:
    #                        ee[3] += 1
    #                        ee[4] = None
    #            # Possible mutation in an end exon
    #            if sum([1 for ee in exons_to_explore if ee[3] < 0]) == len(exons_to_explore):
    #                if not sum([ee[4] for ee in exons_to_explore if not ee[4] is None]):
    #                    #sorted(exons_to_explore, key=lambda x: x[1]-x[2])[0][4] = True
    #                    for ee in exons_to_explore:
    #                        ee[4] = True
    #            # If some are True and some are None get rid of the None
    #            if sum([ee[4] for ee in exons_to_explore if ee[4]]):
    #                if debug:
    #                    sys.stderr.write("getting rid of None " + str(exons_to_explore) + "\n")
    #                exons_to_explore = [ee for ee in exons_to_explore if not ee[4] is None]
    #            # Case where we finished before the end of the exon
    #            if i == len(path)-1:
    #                for ee in exons_to_explore:
    #                    assert ee[4]
    #                    #if not ee[4]:  # Maybe I should make a case for False and None?
    #                    #    log.debug("FINISHING PREMATURELY\n")
    #                    #ee[4] = True
    #            # If all are True
    #            if sum([ee[4] for ee in exons_to_explore if ee[4]]) == len(exons_to_explore):
    #                if debug:
    #                    sys.stderr.write("all true " + str(exons_to_explore) + "\n")
    #                completed = [ee for ee in exons_to_explore if ee[4]]  # necessary?
    #                if debug:
    #                    sys.stderr.write("all true completed " + str(completed) + "\n")
    #                if not exons_in_order:
    #                    # Get all start alternatives from first exon and readjust paths later
    #                    start_orders = set([ee[5] for ee in exons_to_explore])
    #                    if len(exons_to_explore) > 1 and len(start_orders) > 1:
    #                        if debug:
    #                            sys.stderr.write(str(start_orders) + "\n")
    #                        p = []
    #                        for l in range(1, len(path)):
    #                            if path[l] in start_kmers:
    #                                p = path[l:]
    #                                if debug:
    #                                    sys.stderr.write("ALTERNATIVE START\n")
    #                                exon_alternatives.extend(locate_reference(
    #                                    p, exon_seq, exon_names, start_kmers))
    #                                break
    #                choice = sorted(completed, key=lambda x: x[1])[-1]
    #                exons_in_order.append(choice[0])
    #                if debug:
    #                    sys.stderr.write("exons in order " + str(exons_in_order) + "\n")
    #                p = []
    #                rollback = i - abs(choice[2] - choice[1])
    #                for l in range(rollback, len(path)):
    #                    if path[l] in start_kmers:
    #                        if debug:
    #                            sys.stderr.write(str(l) + " " + str(i) + " " + str(rollback) + "\n")
    #                        # Check for overlapping exons
    #                        if l - rollback + 1 <= 31:
    #                            log.debug("Skipping start closer than 31 nt: {}, {}".format(l, rollback))
    #                            continue
    #                            #return None
    #                        p = path[l:]
    #                        break
    #                return find_all(p, exons_in_order, exon_alternatives)
    #        if k in start_kmers:
    #            if debug:
    #                sys.stderr.write("START " + str(k) + "\n")
    #            l = len(exons_to_explore)
    #            for e in exon_num:
    #                if exon_seq[e][0] == k:
    #                    exons_to_explore.append([e, len(exon_seq[e]), 1, 0, False, l])
    #    sys.stderr.write("Loop ended\n")
    #    return None
    #
    #p = find_all(path, [], [])
    #print p
    #return p
