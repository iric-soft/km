#                             -*- Mode: Python -*-
# MutationFinder.py
#

import string
import sys
import logging as log

from . import Graph as ug
from . import PathQuant as upq
from .. utils import common as uc


class MutationFinder:
    def __init__(self, ref_name, ref_seq, jf, graphical, max_stack=500):
        # Load the reference sequence and preparing ref k-mers
        self.first_seq = []
        self.last_seq = []
        self.ref_mer = []
        self.ref_set = []
        for i in range(len(ref_seq)):
            self.first_seq.append(ref_seq[i][0:(jf.k)])
            self.last_seq.append(ref_seq[i][-(jf.k):])
            self.ref_mer.append(uc.get_ref_kmer(ref_seq[i], jf.k, ref_name[i]))
            self.ref_set.append(set(self.ref_mer[i]))
        
        log.debug("Ref. set contains %d kmers.", len(self.ref_set))

        self.ref_seq = ref_seq
        self.jf = jf
        self.node_data = [{} for i in range(len(self.ref_seq))]
        self.done = [set() for i in range(len(self.ref_seq))]
        self.ref_name = ref_name
        
        for i in range(len(ref_seq)):
            self.done[i].add(self.first_seq[i])
            self.node_data[i][self.first_seq[i]] = self.jf.query(self.first_seq[i])
            self.done[i].add(self.last_seq[i])
            self.node_data[i][self.last_seq[i]] = self.jf.query(self.last_seq[i])

        # in case there aren't any
        self.paths = [[] for i in range(len(self.ref_seq))]

        self.max_stack = max_stack

        # register all k-mers from the ref
        for i in range(0, len(self.ref_set)):
            for s in self.ref_set[i]:
                self.node_data[i][s] = self.jf.query(s)

        # kmer walking from each k-mer of ref_seq
        for i in range(len(self.ref_set)):
            self.done[i].update(self.ref_set[i])
            for seq in self.ref_set[i]:
                if seq == self.last_seq:
                    continue
                self.__extend([seq], 0, 0, i)

        
        self.graph_analysis(graphical)

    def __extend(self, stack, breaks, found, target_num):
        """ Recursive depth first search """
        if len(stack) > self.max_stack:
            return
        cur_seq = stack[-1]
        childs = self.jf.get_child(cur_seq, forward=True)

        if len(childs) > 1:
            breaks += 1
            if breaks > 10:
                return

        for child in childs:
            if child in self.done[target_num]:
                self.done[target_num].update(stack)
                self.done[target_num].add(child)
                for p in stack:
                    self.node_data[target_num][p] = self.jf.query(p)
                self.node_data[target_num][cur_seq] = self.jf.query(cur_seq)
                found += 1
            else:
                self.__extend(stack + [child], breaks, found, target_num)

    def graph_analysis(self, graphical=False):
        self.paths = [[] for i in range(len(self.ref_seq))]
        
        kmers = []
        num_k = []
        graph = []
        ref_index = []
        for i in range(len(self.node_data)):
            kmers.append(self.node_data[i].keys())
            num_k.append(len(kmers[i]))
            graph.append(ug.Graph(num_k[i]))
            ref_index.append(list(map(lambda k: kmers[i].index(k), self.ref_mer[i])))

        log.debug("k-mer graph contains %d nodes.", num_k)

        for i in range(len(num_k)):
            for j in range(num_k[i]):
                for k in range(num_k[i]):
                    if j == k:
                        continue
                    if kmers[i][j][1:] == kmers[i][k][:-1]:
                        weight = 1
                        graph[i][j, k] = weight
        
        for l in range(len(ref_index)):
            for k in range(len(ref_index[l])-1):
                i = ref_index[l][k]
                j = ref_index[l][k+1]
                graph[l][i, j] = 0.01

        short_paths = []
        for i in range(len(kmers)):
            graph[i].init_paths(kmers[i].index(self.first_seq[i]),
                            kmers[i].index(self.last_seq[i]))
            short_paths.append(graph[i].all_shortest())

        def get_seq(path, kmer, skip_prefix=True):
            path = list(path)
            if not path:
                # Deals with an empty sequence
                return ""

            if skip_prefix:
                seq = kmer[path[0]][-1]
            else:
                seq = kmer[path[0]]

            for i in path[1:]:
                seq += kmer[i][-1]
            return seq

        def get_name(a, b, target_num, offset=0):
            k = self.jf.k
            diff = graph[target_num].diff_path_without_overlap(a, b, k)
            deletion = diff[3]
            ins = diff[4]

            if (len(a)-len(deletion)+len(ins)) != len(b):
                sys.stderr.write(
                    "ERROR: %s %d != %d" % (
                        "mutation identification could be incorrect",
                        len(a) - len(deletion) + len(ins),
                        len(b)
                    )
                )

                # Fixes cases where we look at two copies of the same sequence
                deletion = diff[3]
                raise Exception()

            # Trim end sequence when in both del and ins:
            del_seq = get_seq(deletion, kmers[target_num], True)
            ins_seq = get_seq(ins, kmers[target_num], True)

            trim = 1
            while (len(del_seq[-trim:]) > 0 and
                    del_seq[-trim:] == ins_seq[-trim:]):
                trim += 1
            trim -= 1
            if trim != 0:
                del_seq = del_seq[:-trim]
                ins_seq = ins_seq[:-trim]

            if diff[0] == diff[1] and not diff[4]:
                return "Reference\t"
            else:
                variant = "Indel"
                # SNP have equal length specific sequences
                if diff[1] == diff[2]:
                    variant = "Substitution"

                # ITD have zero kmers in ref after full trimming.
                # However, this does not distinguish cases where there is
                # garbage between repeats.
                elif diff[0] == diff[5]:
                    variant = "ITD"
                elif len(del_seq) == 0 and len(ins_seq) != 0:
                    variant = "Insertion"
                elif len(del_seq) != 0 and len(ins_seq) == 0:
                    variant = "Deletion"

                return "{}\t{}:{}:{}".format(
                    variant,
                    diff[0] + k + offset,
                    (string.lower(del_seq) + "/" + ins_seq),
                    diff[1] + 1 + offset)

        def get_counts(path, kmer, target_num):
            counts = []
            for i in path:
                counts += [self.node_data[target_num][kmer[i]]]
            # print("length counts: " + str(len(counts)))
            # print("min counts: " + str(min(counts)))

            return counts

        # Quantify all paths independently
        individual = True
        if individual:
            paths_quant = []
            for fusion in range(len(short_paths)):
                for path in short_paths[fusion]:            
                    quant = (upq.PathQuant(all_path=[path, ref_index[fusion]],
                                      counts=self.node_data[fusion].values()))
                    
                    quant.compute_coef()
                    quant.refine_coef()
                    quant.get_ratio()

                    # Reference
                    if list(path) == ref_index[fusion]:
                        quant.adjust_for_reference()

                    paths_quant = quant.get_paths(
                        db_f=self.jf.filename,
                        ref_name= self.ref_name[fusion],
                        name_f=lambda path: get_name(ref_index[fusion], path, fusion),
                        seq_f=lambda path: get_seq(path, kmers[fusion], skip_prefix=False),
                        ref_path=ref_index[fusion], info="vs_ref",
                        get_min_f=lambda path: min(get_counts(path, kmers[fusion], fusion)))
                    
                    
                    self.paths[fusion] += paths_quant


                if graphical:
                    import matplotlib.pyplot as plt

                    plt.figure(figsize=(10, 6))
                    for path in short_paths:
                            plt.plot(get_counts(path, kmers[fusion], fusion),
                                label=get_name(ref_index[fusion], path).replace("\t", " "))
                    plt.legend()
                    plt.show()

        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
    '''    cluster = True
        if cluster:
            variant_diffs = []
            diff = []
            variant_diffs = []
            for fusion in range(len(short_paths)):
                variant_set = set(range(0, len(short_paths[fusion])))
                for variant in short_paths[fusion]:
                    diff.append(graph[fusion].diff_path_without_overlap(
                        ref_index[fusion], variant, self.jf.k))
                    variant_diffs += [diff[fusion]]

                def get_intersect(start, stop):
                    for var in variant_set:
                        if (variant_diffs[var][1] >= start and
                                variant_diffs[var][0] <= stop):
                            return var
                    return -1

                variant_groups = []
                while len(variant_set) > 0:
                    seed = variant_set.pop()
                    grp = [seed]
                    start = variant_diffs[seed][0]
                    stop = variant_diffs[seed][1]
                    variant = get_intersect(start, stop)
                    while variant != -1:
                        variant_set.remove(variant)
                        grp += [variant]
                        start = min(start, variant_diffs[variant][0])
                        stop = max(stop, variant_diffs[variant][1])
                        variant = get_intersect(start, stop)
                    variant_groups += [(start, stop, grp)]

                num_cluster = 0
                for var_gr in variant_groups:
                    if (len(var_gr[2]) == 1 and
                            list(short_paths[var_gr[2][0]]) == ref_index):
                        continue
                    num_cluster += 1

                    start = var_gr[0]
                    stop = var_gr[1]
                    var_size = max([abs(x[2]-x[1]+1) for x in [variant_diffs[v] for v in var_gr[2]]])
                    offset = max(0, start - var_size)
                    ref_path = ref_index[offset:stop]
                    clipped_paths = [ref_path]
                    for var in var_gr[2]:
                        start_off = offset
                        stop_off = variant_diffs[var][2] + (stop - variant_diffs[var][1])
                        clipped_paths += [short_paths[var][start_off:stop_off]]

                    quant = upq.PathQuant(all_path=clipped_paths,
                                      counts=self.node_data.values())

                    quant.compute_coef()
                    quant.refine_coef()

                    quant.get_ratio()

                    for i in range(self.ref_name):
                        paths_quant.append(quant.get_paths(
                        db_f=self.jf.filename,
                        ref_name=self.ref_name[i],
                        name_f=lambda path: get_name(ref_path, path, offset),
                        seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                        ref_path=ref_path,
                        info="cluster %d n=%d" % (num_cluster, len(var_gr[2])),
                        get_min_f=lambda path: min(get_counts(path, kmer, i)),
                        start_off=start_off))
                        paths_quant[i]
                        self.paths[i] += paths_quant[i]

                    if graphical:
                        import matplotlib.pyplot as plt

                        plt.figure(figsize=(10, 6))
                        for i in range(self.get_paths):
                            for path, ratio in zip(clipped_paths, quant.get_ratio()):
                                if path == ref_path:
                                    plt.plot(get_counts(path, kmer, i),
                                        label="Reference")
                                else:
                                    plt.plot(get_counts(path, kmer, i),
                                        label=get_name(ref_path, path, offset).split("\t")[0])
                        plt.legend()
                        plt.show()
'''
    def get_paths(self, target_num):
        return self.paths[target_num]

    @staticmethod
    def output_header():
        upq.PathQuant.output_header()
