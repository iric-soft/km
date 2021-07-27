#                             -*- Mode: Python -*-
# MutationFinder.py
#

import sys
import logging as log

from . import Graph as ug
from . import PathQuant as upq
from .. utils import common as uc


class MutationFinder:
    def __init__(self, ref_name, ref_seq, jf, graphical, max_stack=500,
                 max_break=10):
        # Load the reference sequence and preparing ref k-mers

        self.first_seq = ref_seq[0:(jf.k)]
        self.last_seq = ref_seq[-(jf.k):]

        self.ref_mer = uc.get_ref_kmer(ref_seq, jf.k, ref_name)
        self.ref_set = set(self.ref_mer)
        log.debug("Ref. set contains %d kmers.", len(self.ref_set))

        self.ref_seq = ref_seq
        self.jf = jf
        self.node_data = {}
        self.done = set()
        self.ref_name = ref_name

        self.done.add(self.first_seq)
        self.node_data[self.first_seq] = self.jf.query(self.first_seq)
        self.done.add(self.last_seq)
        self.node_data[self.last_seq] = self.jf.query(self.last_seq)

        # in case there aren't any
        self.paths = []

        self.max_stack = max_stack
        self.max_break = max_break

        # register all k-mers from the ref
        for s in self.ref_set:
            self.node_data[s] = self.jf.query(s)

        # kmer walking from each k-mer of ref_seq
        self.done.update(self.ref_set)
        for seq in self.ref_set:
            if seq == self.last_seq:
                continue
            self.__extend([seq], 0, 0)

        self.kmer = list(self.node_data.keys())
        self.kmer_count = list(self.node_data.values())
        self.num_k = len(self.kmer)

        # reference path, with node indices
        self.ref_index = [self.kmer.index(k) for k in self.ref_mer]

        log.debug("k-mer graph contains %d nodes.", self.num_k)

        self.graph_analysis(graphical)

    def __extend(self, stack, breaks, found):
        """ Recursive depth first search """
        if len(stack) > self.max_stack:
            return
        cur_seq = stack[-1]
        childs = self.jf.get_child(cur_seq, forward=True)

        if len(childs) > 1:
            breaks += 1
            if breaks > self.max_break:
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

    def _diff_path_without_overlap(self, ref, seq):
        # Returns (start, stop_ref, stop_variant, kmers_ref, kmers_variant, stop_ref_fully_trimmed)
        i = 0

        while i < len(ref) and i < len(seq) and ref[i] == seq[i]:
            i += 1

        j_ref = len(ref)
        j_seq = len(seq)
        while j_ref > i + (self.jf.k - 1) and j_seq > i + (self.jf.k - 1) and ref[j_ref - 1] == seq[j_seq - 1]: #  + (k - 1):to prevent kmer from overlapping
            j_ref -= 1
            j_seq -= 1

        k_ref = j_ref
        k_seq = j_seq
        while k_ref > i and ref[k_ref - 1] == seq[k_seq - 1]:
            k_ref -= 1
            k_seq -= 1

        # log.debug("diffpath : " + " ".join (str(x) for x in [i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref]))

        return (i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref)

    def get_seq(self, path, skip_prefix=True):
        path = list(path)

        if not path:
            # Deals with an empty sequence
            return ""

        if skip_prefix:
            seq = self.kmer[path[0]][-1]
        else:
            seq = self.kmer[path[0]]

        for i in path[1:]:
            seq += self.kmer[i][-1]

        return seq

    def get_name(self, ref_ix, path_ix, offset=0):
        diff = self._diff_path_without_overlap(ref_ix, path_ix)
        deletion = diff[3]
        ins = diff[4]

        if (len(ref_ix)-len(deletion)+len(ins)) != len(path_ix):
            sys.stderr.write(
                "ERROR: %s %d != %d" % (
                    "mutation identification could be incorrect",
                    len(ref_ix) - len(deletion) + len(ins),
                    len(path_ix)
                )
            )

            # Fixes cases where we look at two copies of the same sequence
            deletion = diff[3]
            raise Exception()

        # Trim end sequence when in both del and ins:
        del_seq = self.get_seq(deletion, True)
        ins_seq = self.get_seq(ins, True)

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
                diff[0] + self.jf.k + offset,
                (str.lower(del_seq) + "/" + ins_seq),
                diff[1] + 1 + offset)

    def get_counts(self, path):
        counts = []
        for i in path:
            counts += [self.node_data[self.kmer[i]]]
        # print("length counts: " + str(len(counts)))
        # print("min counts: " + str(min(counts)))

        return counts

    def graph_analysis(self, graphical=False):
        self.paths = []
        kmer = self.kmer

        num_k = self.num_k
        graph = ug.Graph(num_k)
        # The reference path, with node numbers
        ref_index = self.ref_index

        for i in range(num_k):
            for j in range(num_k):
                if i == j:
                    continue
                if kmer[i][1:] == kmer[j][:-1]:
                    weight = 1
                    graph[i, j] = weight

        for k in range(len(ref_index)-1):
            i = ref_index[k]
            j = ref_index[k+1]
            graph[i, j] = 0.01

        graph.init_paths(kmer.index(self.first_seq),
                         kmer.index(self.last_seq))

        self.short_paths = graph.all_shortest()

        individual = True
        self.quantify_paths(graphical, do=individual)

        cluster = True
        self.find_clusters(graphical, do=cluster)

    def quantify_paths(self, graphical=False, do=True):
        # Quantify all paths independently
        if do:
            ref_index = self.ref_index
            for path in self.short_paths:
                quant = upq.PathQuant(all_paths=[path, ref_index],
                                      counts=list(self.node_data.values()))

                quant.compute_coef()
                quant.refine_coef()
                quant.get_ratio()

                # Reference
                if list(path) == ref_index:
                    quant.adjust_for_reference()

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: self.get_name(ref_index, path),
                    seq_f=lambda path: self.get_seq(path, skip_prefix=False),
                    ref_path=ref_index, info="vs_ref",
                    get_min_f=lambda path: min(self.get_counts(path)))

                self.paths += self.paths_quant

            if graphical:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(10, 6))
                for path in self.short_paths:
                    plt.plot(self.get_counts(path),
                             label=self.get_name(ref_index, path).replace("\t", " "))
                plt.legend()
                plt.show()

    def find_clusters(self, graphical=False, do=True):
        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        if do:
            ref_index = self.ref_index
            variant_diffs = []
            variant_set = set(range(0, len(self.short_paths)))
            for variant in self.short_paths:
                diff = self._diff_path_without_overlap(
                    self.ref_index, variant)
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
                        list(self.short_paths[var_gr[2][0]]) == ref_index):
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
                    clipped_paths += [self.short_paths[var][start_off:stop_off]]

                quant = upq.PathQuant(all_paths=clipped_paths,
                                      counts=list(self.node_data.values()))

                quant.compute_coef()
                quant.refine_coef()

                quant.get_ratio()

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: self.get_name(ref_path, path, offset),
                    seq_f=lambda path: self.get_seq(path, skip_prefix=False),
                    ref_path=ref_path,
                    info="cluster %d n=%d" % (num_cluster, len(var_gr[2])),
                    get_min_f=lambda path: min(self.get_counts(path)),
                    start_off=start_off)

                self.paths_quant

                self.paths += self.paths_quant

                if graphical:
                    import matplotlib.pyplot as plt

                    plt.figure(figsize=(10, 6))
                    for path, ratio in zip(clipped_paths, quant.get_ratio()):
                        if path == ref_path:
                            plt.plot(self.get_counts(path),
                                     label="Reference")
                        else:
                            plt.plot(self.get_counts(path),
                                     label=self.get_name(ref_path, path, offset).split("\t")[0])
                    plt.legend()
                    plt.show()

    def get_paths(self, sort=True):
        if sort:
            self.paths = sorted(self.paths, key=lambda x: (x[3], x[2], x[6]))
        return self.paths

    def get_paths_quant(self):
        return self.paths_quant

    @staticmethod
    def output_header():
        upq.Path.output_header()
