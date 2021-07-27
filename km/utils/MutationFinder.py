#                             -*- Mode: Python -*-
# MutationFinder.py
#

import sys
import logging as log
from collections import namedtuple

from . import Graph as ug
from . import PathQuant as upq
from .. utils import common as uc


PathDiff = namedtuple('PathDiff',
    [
        'start',
        'end_ref',
        'end_var',
        'kmers_ref',
        'kmers_var',
        'end_ref_overlap'
    ])


class MutationFinder:
    def __init__(self, ref_name, ref_seq, jf, max_stack=500,
                 max_break=10):
        # Load the reference sequence and preparing ref k-mers

        self.first_seq = ref_seq[0:(jf.k)]
        self.last_seq = ref_seq[-(jf.k):]

        self.ref_mer = uc.get_ref_kmer(ref_seq, jf.k, ref_name)
        self.ref_set = set(self.ref_mer)
        log.debug("Ref. set contains %d kmers.", len(self.ref_set))

        self.ref_seq = ref_seq
        self.ref_name = ref_name
        self.jf = jf
        self.max_stack = max_stack
        self.max_break = max_break

        self.node_data = {}
        self.done = set()

        self.done.add(self.first_seq)
        self.node_data[self.first_seq] = self.jf.query(self.first_seq)
        self.done.add(self.last_seq)
        self.node_data[self.last_seq] = self.jf.query(self.last_seq)

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

    def _diff_path_without_overlap(self, ref, seq, k):
        # Returns (start, stop_ref, stop_variant, kmers_ref, kmers_variant, stop_ref_fully_trimmed)

        i = 0
        while True:
            if i == len(ref):
                break
            if i == len(seq):
                break
            if ref[i] != seq[i]:
                break
            # all conditions pass, can increase position by 1
            i += 1
        assert i <= len(ref) and i <= len(seq)

        j_ref = len(ref)
        j_seq = len(seq)
        # Note:
        #   - `+ k` == to prevent kmers from overlapping
        while True:
            if j_ref < i + k:
                break
            if j_seq < i + k:
                break
            if ref[j_ref - 1] != seq[j_seq - 1]:
                break
            # all conditions pass, can decrease positions by 1
            j_ref -= 1
            j_seq -= 1

        k_ref = j_ref
        k_seq = j_seq
        while True:
            if k_ref <= i:
                break
            if ref[k_ref - 1] != seq[k_seq - 1]:
                break
            # all conditions pass, can decrease positions by 1
            k_ref -= 1
            k_seq -= 1

        path_diff = PathDiff(i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref)

        return path_diff

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
        diff = self._diff_path_without_overlap(ref_ix, path_ix, self.jf.k)

        if len(ref_ix) - len(diff.kmers_ref) + len(diff.kmers_var) != len(path_ix):
            sys.stderr.write(
                "ERROR: %s %d != %d" % (
                    "mutation identification could be incorrect",
                    len(ref_ix) - len(diff.kmers_ref) + len(diff.kmers_var),
                    len(path_ix)
                )
            )
            # Fixes cases where we look at two copies of the same sequence
            raise Exception()

        # Get sequence from indices
        del_seq = self.get_seq(diff.kmers_ref, skip_prefix=True)
        ins_seq = self.get_seq(diff.kmers_var, skip_prefix=True)

        trim = 1
        if len(del_seq) > 0:
            assert del_seq != ins_seq
            while del_seq[-trim:] == ins_seq[-trim:]:
                trim += 1
        trim -= 1

        if trim != 0:
            del_seq = del_seq[:-trim]
            ins_seq = ins_seq[:-trim]

        if diff.end_ref == diff.end_var:
            if diff.start == diff.end_ref:
                # Ref and Alt end on at the same position, which is also the start
                variant = "Reference"
            else:
                # SNPs and MNPs come from sequences of equal lengths
                variant = "Substitution"
        else:
            if diff.start == diff.end_ref_overlap:
                # we retraced the whole reference
                # ITD have zero kmers in ref after trimming, but I&I do.
                variant = "ITD"
            else:
                # we have an indel
                variant = "Indel"
                if diff.end_ref < diff.end_var:
                    if len(del_seq) == 0:
                        variant = "Insertion"
                else:
                    if len(ins_seq) == 0:
                        variant = "Deletion"

        if variant == "Reference":
            return variant + "\t"
        else:
            return "{}\t{}:{}:{}".format(
                variant,
                diff.start + self.jf.k + offset,
                str.lower(del_seq) + "/" + ins_seq,
                diff.end_ref + 1 + offset
            )

    def get_counts(self, path):
        counts = [self.node_data[self.kmer[k]] for k in path]
        return counts

    def graph_analysis(self, graphical=False):
        self.paths = []

        # Initialize graph
        graph = ug.Graph(self.num_k)

        # For all kmers, find the next kmers with k - 1 overlap
        # and assign a weight of 1
        weight = 1
        for i in range(self.num_k):
            for j in range(self.num_k):
                if i != j:
                    if self.kmer[i][1:] == self.kmer[j][:-1]:
                        graph[i, j] = weight

        # Change the weights for reference from 1 to 0.01 in the graph
        weight = 0.01
        for k in range(len(self.ref_index)-1):
            i = self.ref_index[k]
            j = self.ref_index[k+1]
            graph[i, j] = weight

        # Initialize paths and remove reference edges
        graph.init_paths(
            self.kmer.index(self.first_seq),
            self.kmer.index(self.last_seq)
        )

        self.short_paths = graph.all_shortest()

    def quantify_paths(self, graphical=False):
        # Quantify all paths independently
        if graphical:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(10, 6))
            for path in self.short_paths:
                plt.plot(
                    self.get_counts(path),
                    label=self.get_name(
                            self.ref_index,
                            path
                        ).replace("\t", " ")
                )
            plt.legend()
            plt.show()

        for path in self.short_paths:
            quant = upq.PathQuant(
                all_paths=[path, self.ref_index],
                counts=self.kmer_count
            )
            quant.compute_coef()
            quant.refine_coef()
            quant.get_ratio()

            # Reference
            if list(path) == self.ref_index:
                quant.adjust_for_reference()

            rvaf, ref_rvaf = quant.rVAF
            coef, ref_coef = quant.coef

            path_o = upq.Path(
                self.jf.filename,
                self.ref_name,
                self.get_name(self.ref_index, path),
                rvaf,
                coef,
                min(self.get_counts(path)),
                0,
                self.get_seq(path, skip_prefix=False),
                ref_rvaf,
                ref_coef,
                self.get_seq(self.ref_index, skip_prefix=False),
                "vs_ref"
            )

            self.paths.append(path_o)

    def find_clusters(self, graphical=False):
        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        variant_diffs = []
        variant_set = set(range(0, len(self.short_paths)))
        for variant in self.short_paths:
            diff = self._diff_path_without_overlap(self.ref_index, variant, self.jf.k)
            variant_diffs.append(diff)

        def get_intersect(start, stop):
            for var in variant_set:
                cur_start = variant_diffs[var].start
                cur_end = variant_diffs[var].end_ref
                if cur_end >= start and cur_start <= stop:
                    return var
            return -1

        variant_groups = []
        while len(variant_set) > 0:
            seed = variant_set.pop()
            grp = [seed]
            start = variant_diffs[seed].start
            stop = variant_diffs[seed].end_ref
            variant = get_intersect(start, stop)
            while variant != -1:
                variant_set.remove(variant)
                grp += [variant]
                start = min(start, variant_diffs[variant].start)
                stop = max(stop, variant_diffs[variant].end_ref)
                variant = get_intersect(start, stop)
            variant_groups.append((start, stop, grp))

        self.clusters = []

        for var_gr in variant_groups:
            start, stop, grp_ixs = var_gr

            if len(grp_ixs) == 1:
                var = grp_ixs[0]
                path_index = list(self.short_paths[var])
                if path_index == self.ref_index:
                    continue

            var_diffs = [variant_diffs[v] for v in grp_ixs]
            var_size = max(
                [abs(d.end_var - d.end_ref + 1) for d in var_diffs]
            )
            offset = max(0, start - var_size)
            ref_path = self.ref_index[offset:stop]

            clipped_paths = []
            for var in grp_ixs:
                cur_diff = variant_diffs[var]
                stop_off = cur_diff.end_var + stop - cur_diff.end_ref
                new_path = self.short_paths[var][offset:stop_off]
                clipped_paths.append(new_path)

            self.clusters.append((ref_path, clipped_paths, offset))

        if graphical:
            import matplotlib.pyplot as plt

        for i, cluster in enumerate(self.clusters):
            ref_path, clipped_paths, start_off = cluster
            num_cluster = i + 1

            if graphical:
                plt.figure(figsize=(10, 6))
                plt.plot(
                    self.get_counts(ref_path),
                    label="Reference"
                )
                for path in clipped_paths:
                    assert path != ref_path
                    plt.plot(
                        self.get_counts(path),
                        label=self.get_name(
                            ref_path,
                            path,
                            start_off
                        ).split("\t")[0]
                    )
                plt.legend()
                plt.show()

            quant = upq.PathQuant(
                all_paths=[ref_path] + clipped_paths,
                counts=self.kmer_count
            )
            quant.compute_coef()
            quant.refine_coef()
            quant.get_ratio()

            ref_rvaf, paths_rvaf = quant.rVAF[0], quant.rVAF[1:]
            ref_coef, paths_coef = quant.coef[0], quant.coef[1:]

            for path, rvaf, coef in zip(clipped_paths, paths_rvaf, paths_coef):
                path_o = upq.Path(
                    self.jf.filename,
                    self.ref_name,
                    self.get_name(ref_path, path, start_off),
                    rvaf,
                    coef,
                    min(self.get_counts(path)),
                    start_off,
                    self.get_seq(path, skip_prefix=False),
                    ref_rvaf,
                    ref_coef,
                    self.get_seq(ref_path, skip_prefix=False),
                    "cluster %d n=%d" % (num_cluster, len(clipped_paths))
                )

                self.paths.append(path_o)

    def get_paths(self, sort=True):
        if sort:
            paths = sorted(self.paths, key=lambda x: (x[3], x[2], x[6]))
        else:
            paths = self.paths

        return paths

    @staticmethod
    def output_header():
        upq.Path.output_header()
