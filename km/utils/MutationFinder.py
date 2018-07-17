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
    def __init__(self, target_name, target_seq, jf, graphical, max_stack=500,
                 max_break=10):
        # Load the reference sequence and preparing ref k-mers
        self.first_kmer = []
        self.last_kmer = []
        self.target_kmers = []
        self.target_set = []
        sum_kmers = 0
        for i in range(len(target_seq)):
            self.first_kmer.append(target_seq[i][0:(jf.k)])
            self.last_kmer.append(target_seq[i][-(jf.k):])
            self.target_kmers.append(uc.get_target_kmers(target_seq[i], jf.k, target_name[i]))
            self.target_set.append(set(self.target_kmers[i]))
            sum_kmers += len(self.target_kmers[i])

        log.debug("Ref. set contains %d kmers.", sum_kmers)

        self.target_seq = target_seq
        self.jf = jf
        self.node_data = {}
        self.done = [set() for i in range(len(target_seq))]
        self.target_name = target_name

        for i in range(len(self.first_kmer)):
            self.done[i].add(self.first_kmer[i])
            self.node_data[self.first_kmer[i]] = self.jf.query(self.first_kmer[i])
            self.done[i].add(self.last_kmer[i])
            self.node_data[self.last_kmer[i]] = self.jf.query(self.last_kmer[i])

        # in case there aren't any
        self.paths = []

        self.max_stack = max_stack
        self.max_break = max_break

        # register all k-mers from the ref
        for target in range(len(self.target_set)):
            # kmer walking from each k-mer of target_seq
            self.done[target].update(self.target_set[target])
            for kmer in self.target_set[target]:
                self.node_data[kmer] = self.jf.query(kmer)

        for i in range(len(self.last_kmer)):
            for kmer in self.target_set[i]:
                if kmer == self.last_kmer[i]:
                    continue
                self.__extend([kmer], 0, i)

        self.graph_analysis(graphical)

    def __extend(self, stack, breaks, target):
        """ Recursive depth first search """
        if len(stack) > self.max_stack:
            return
        cur_kmer = stack[-1]
        childs = self.jf.get_child(cur_kmer, forward=True)

        if len(childs) > 1:
            breaks += 1
            if breaks > self.max_break:
                return

        for child in childs:
            if child in self.done[target]:
                self.done[target].update(stack)
                for p in stack:
                    self.node_data[p] = self.jf.query(p)
            else:
                self.__extend(stack + [child], breaks, target)

    def graph_analysis(self, graphical=False):
        self.paths = []
        kmer = self.node_data.keys()

        num_k = len(kmer)
        graph = ug.Graph(num_k)
        # The reference path, with node numbers
        target_index = []
        for i in range(len(self.target_kmers)):
            target_index.append(map(lambda k: kmer.index(k), self.target_kmers[i]))

        log.debug("k-mer graph contains %d nodes.", num_k)

        for i in range(num_k):
            for j in range(num_k):
                if i == j:
                    continue
                if kmer[i][1:] == kmer[j][:-1]:
                    weight = 1
                    graph[i, j] = weight

        for l in range(len(target_index)):
            for k in range(len(target_index[l])-1):
                i = target_index[l][k]
                j = target_index[l][k+1]
                graph[i, j] = 0.01

        short_paths = []
        for i in range(len(self.first_kmer)):
            graph.init_paths(kmer.index(self.first_kmer[i]),
                             kmer.index(self.last_kmer[i]))

            short_paths.append(graph.all_shortest())
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

        def get_name(a, b, offset=0):
            k = self.jf.k
            diff = graph.diff_path_without_overlap(a, b, k)
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
            del_seq = get_seq(deletion, kmer, True)
            ins_seq = get_seq(ins, kmer, True)

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

        def get_counts(path, kmer):
            counts = []
            for i in path:
                counts += [self.node_data[kmer[i]]]
            # print("length counts: " + str(len(counts)))
            # print("min counts: " + str(min(counts)))

            return counts

        # Quantify all paths independently
        individual = True
        if individual:
            for target_id in range(len(short_paths)):
                for path in short_paths[target_id]:
                    print path
                    quant = upq.PathQuant(all_path=[path, target_index[target_id]],
                                          counts=self.node_data.values())

                    quant.compute_coef()
                    quant.refine_coef()
                    quant.get_ratio()

                    # Reference
                    if list(path) == target_index[target_id]:
                        quant.adjust_for_reference()

                    paths_quant = quant.get_paths(
                        db_f=self.jf.filename,
                        target_name=self.target_name[target_id],
                        name_f=lambda path: get_name(target_index[target_id], path),
                        seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                        target_path=target_index[target_id], info="vs_ref",
                        get_min_f=lambda path: min(get_counts(path, kmer)))

                    self.paths += paths_quant

                if graphical:
                    import matplotlib.pyplot as plt

                    plt.figure(figsize=(10, 6))
                    for path in short_paths:
                        plt.plot(get_counts(path, kmer),
                                 label=get_name(target_index[target_id], path).replace("\t", " "))
                    plt.legend()
                    plt.show()

        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        #
        cluster = True
        if cluster:
            for target_id in range(len(short_paths)):
                variant_set = set(range(0, len(short_paths[target_id])))
                variant_diffs = []
                for variant in short_paths[target_id]:
                    diff = graph.diff_path_without_overlap(
                        target_index[target_id], variant, self.jf.k)
                    variant_diffs += [diff]

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
                          list(short_paths[target_id][var_gr[2][0]]) == target_index[target_id]):
                        continue
                    num_cluster += 1

                    start = var_gr[0]
                    stop = var_gr[1]
                    var_size = max([abs(x[2]-x[1]+1) for x in [variant_diffs[v] for v in var_gr[2]]])
                    offset = max(0, start - var_size)
                    target_path = target_index[target_id][offset:stop]
                    clipped_paths = [target_path]
                    for var in var_gr[2]:
                        start_off = offset
                        stop_off = variant_diffs[var][2] + (stop - variant_diffs[var][1])
                        clipped_paths += [short_paths[target_id][var][start_off:stop_off]]

                    quant = upq.PathQuant(all_path=clipped_paths,
                                      counts=self.node_data.values())

                    quant.compute_coef()
                    quant.refine_coef()

                    quant.get_ratio()

                    paths_quant = (quant.get_paths(
                            db_f=self.jf.filename,
                            target_name=self.target_name[target_id],
                            name_f=lambda path: get_name(target_path, path, offset),
                            seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                            target_path=target_path,
                            info="cluster %d n=%d" % (num_cluster, len(var_gr[2])),
                            get_min_f=lambda path: min(get_counts(path, kmer )),
                            start_off=start_off))

                    paths_quant
                    self.paths += paths_quant

                    if graphical:
                        import matplotlib.pyplot as plt

                        plt.figure(figsize=(10, 6))
                        for i in range(self.get_paths):
                            for path, ratio in zip(clipped_paths, quant.get_ratio()):
                                if path == target_path:
                                    plt.plot(get_counts(path, kmer ),
                                        label="Reference")
                                else:
                                    plt.plot(get_counts(path, kmer ),
                                        label=get_name(target_path, path, offset).split("\t")[0])
                        plt.legend()
                        plt.show()

    def get_paths(self):
        return self.paths

    @staticmethod
    def output_header():
        upq.PathQuant.output_header()
