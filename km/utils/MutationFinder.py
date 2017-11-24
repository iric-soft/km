#                             -*- Mode: Python -*-
# MutationFinder.py
#

import string
import sys
import logging as log

from . import Graph as ug
from . import PathQuant as upq


class MutationFinder:
    def __init__(self, ref_name, ref_seq, jf, graphical, max_stack=500):
        # Load the reference sequence and preparing ref k-mers

        self.first_seq = ref_seq[0:(jf.k)]
        self.last_seq = ref_seq[-(jf.k):]

        # Load reference kmers.
        ref_mer = []
        for i in range(len(ref_seq) - (jf.k) + 1):
            kmer = ref_seq[i:(i + jf.k)]
            if kmer in ref_mer:
                sys.stderr.write(
                    "WARNING: %s found multiple times in reference %s %d" % (
                        kmer, ref_name, i
                    )
                )

            ref_mer.append(kmer)
        ref_set = set(ref_mer)

        log.debug("Ref. set contains %d kmers.", len(ref_set))

        self.ref_set = ref_set
        self.ref_seq = ref_seq
        self.jf = jf
        self.node_data = {}
        self.done = set()
        self.ref_mer = ref_mer
        self.ref_name = ref_name

        self.done.add(self.first_seq)
        self.node_data[self.first_seq] = self.jf.query(self.first_seq)
        self.done.add(self.last_seq)
        self.node_data[self.last_seq] = self.jf.query(self.last_seq)

        # in case there aren't any
        self.paths = []

        self.max_stack = max_stack

        # register all k-mers from the ref
        for s in self.ref_set:
            self.node_data[s] = self.jf.query(s)

        # kmer walking from each k-mer of ref_seq
        self.done.update(self.ref_set)
        for seq in ref_set:
            if seq == self.last_seq:
                continue
            self.__extend([seq], 0, 0)

        self.graph_analysis(graphical)

    def __extend(self, stack, breaks, found):
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
            if child in self.done:
                self.done.update(stack)
                self.done.add(child)
                for p in stack:
                    self.node_data[p] = self.jf.query(p)
                self.node_data[cur_seq] = self.jf.query(cur_seq)
                found += 1
            else:
                self.__extend(stack + [child], breaks, found)

    def graph_analysis(self, graphical=False):
        self.paths = []
        kmer = self.node_data.keys()

        k_len = len(kmer)
        graph = ug.Graph(k_len)
        # The reference path, with node numbers
        ref_i = map(lambda k: kmer.index(k), self.ref_mer)

        log.debug("k-mer graph contains %d nodes.", k_len)

        for i in range(k_len):
            for j in range(k_len):
                if i == j:
                    continue
                if kmer[i][1:] == kmer[j][:-1]:
                    weight = 1
                    graph[i, j] = weight

        for k in range(len(ref_i)-1):
            i = ref_i[k]
            j = ref_i[k+1]
            graph[i, j] = 0.01

        graph.init_paths(kmer.index(self.first_seq),
                         kmer.index(self.last_seq))
        short_paths = graph.all_shortest()

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
                    "WARNING: %s %d != %d" % (
                        "mutation identification could be incorrect",
                        len(a) - len(deletion) + len(ins),
                        len(b)
                    )
                )

                # Fixes cases where we look at two copies of the same sequence
                deletion = diff[3]
                raise Exception()

            old_del = diff[3][:-(k - 1)]
            old_ins = diff[4][:-(k - 1)]
            old_del_seq = get_seq(old_del, kmer, True)
            old_ins_seq = get_seq(old_ins, kmer, True)

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

            if old_del_seq != del_seq or old_ins_seq != ins_seq:
                sys.stderr.write("### WARNING: INCOHERENCE IN NAMING ###\n")
                sys.stderr.write("DEL %s %s\n" % (old_del_seq, del_seq))
                sys.stderr.write("INS %s %s\n" % (old_ins_seq, ins_seq))
                sys.stderr.write("### END WARNING ###\n")
                # raise Exception()

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
            print("length counts: " + str(len(counts)))
            print("min counts: " + str(min(counts)))

            return counts

        # Consider all paths at once
        whole = False
        if whole:
            quant = upq.PathQuant(all_path=short_paths,
                                  counts=self.node_data.values())

            quant.compute_coef()
            quant.refine_coef()
            quant.get_ratio()

            self.paths_quant = quant.get_paths(
                db_f=self.jf.filename,
                ref_name=self.ref_name,
                name_f=lambda path: get_name(ref_i, path),
                seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                ref_path=ref_i, info="whole sequence",
                get_min_f=lambda path: min(get_counts(path, kmer)))

            self.paths += self.paths_quant

            # if graphical:
            #     import matplotlib.pyplot as plt

            #     plt.figure(figsize=(10, 6))
            #     for path, ratio in zip(short_paths, quant.get_ratio()):
            #         plt.plot(get_counts(path, kmer),# * ratio,
            #                   label=get_name(ref_i, path).split("\t")[0])
            #     plt.legend()
            #     plt.show()

        # Quantify all paths independently
        individual = True
        if individual:
            for path in short_paths:
                quant = upq.PathQuant(all_path=[path, ref_i],
                                      counts=self.node_data.values())

                quant.compute_coef()
                quant.refine_coef()
                quant.get_ratio()

                # Reference
                if list(path) == ref_i:
                    quant.adjust_for_reference()

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name(ref_i, path),
                    seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                    ref_path=ref_i, info="vs_ref",
                    get_min_f=lambda path: min(get_counts(path, kmer)))

                self.paths += self.paths_quant

            if graphical:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(10, 6))
                for path in short_paths:
                    plt.plot(get_counts(path, kmer),
                             label=get_name(ref_i, path).replace("\t", " "))
                plt.legend()
                plt.show()

        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        cluster = True
        if cluster:
            if len(short_paths) == 1 and list(short_paths[0]) == ref_i:
                quant = upq.PathQuant(all_path=short_paths,
                                      counts=self.node_data.values())

                quant.compute_coef()
                quant.refine_coef()

                quant.get_ratio()

                if 0 in get_counts(short_paths[0], kmer):
                    quant.ratio = quant.coef * 0

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name(ref_i, path),
                    seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                    ref_path=short_paths, info="cluster 0",
                    get_min_f=lambda path: min(get_counts(path, kmer)))

                self.paths += self.paths_quant
            else:
                variant_diffs = []
                variant_set = set(range(0, len(short_paths)))
                for variant in short_paths:
                    diff = graph.diff_path_without_overlap(
                        ref_i, variant, self.jf.k)
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
                            list(short_paths[var_gr[2][0]]) == ref_i):
                        continue
                    num_cluster += 1

                    start = var_gr[0]
                    stop = var_gr[1]
                    var_size = max([abs(x[2]-x[1]+1) for x in [variant_diffs[v] for v in var_gr[2]]])
                    offset = max(0, start - var_size)
                    ref_path = ref_i[offset:stop]
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

                    self.paths_quant = quant.get_paths(
                        db_f=self.jf.filename,
                        ref_name=self.ref_name,
                        name_f=lambda path: get_name(ref_path, path, offset),
                        seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                        ref_path=ref_path,
                        info="cluster %d n=%d" % (num_cluster, len(var_gr[2])),
                        get_min_f=lambda path: min(get_counts(path, kmer)))

                    self.paths += self.paths_quant

                    if graphical:
                        import matplotlib.pyplot as plt

                        plt.figure(figsize=(10, 6))
                        for path, ratio in zip(clipped_paths, quant.get_ratio()):
                            if path == ref_path:
                                plt.plot(get_counts(path, kmer),
                                         label="Reference")
                            else:
                                plt.plot(get_counts(path, kmer),
                                         label=get_name(ref_path, path, offset).split("\t")[0])
                        plt.legend()
                        plt.show()

    @staticmethod
    def output_header():
        upq.PathQuant.output_header()

    def get_paths(self):
        return self.paths

    def get_paths_quant(self):
        return self.paths_quant
