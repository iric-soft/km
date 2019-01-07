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
    def __init__(self, ref_name, ref_seq, jf, graphical, max_stack=500,
                 max_break=10, mode='mutation'):
        
        # Load the reference sequence and prepare ref k-mers
        self.first_kmer = "BigBang"
        self.last_kmer = "BigCrunch"
        
        if mode == "fusion":
            self.start_kmers = set([seq[0:(jf.k)] for seq in ref_seq[0]])
            self.end_kmers = set([seq[-(jf.k):] for seq in ref_seq[1]])
            
            ref_mers_left = [[self.first_kmer] + uc.get_ref_kmers(seq, jf.k, name)
                             for seq, name in zip(ref_seq[0], ref_name[0])]
            ref_mers_right = [uc.get_ref_kmers(seq, jf.k, name) + [self.last_kmer]
                              for seq, name in zip(ref_seq[1], ref_name[1])]
            self.ref_mers = ref_mers_left + ref_mers_right
            self.ref_set = set([kmer for seq in self.ref_mers for kmer in seq])
            
            multiple = set([k for seq in ref_mers_left for k in seq]).intersection(
                    set([k for seq in ref_mers_right for k in seq]))
            if multiple:
                raise ValueError("%s found on both extremities" % (", ".join(multiple)))
            
            self.ref_name_fusion = ref_name
            ref_name = ref_name[0][0].split("_")[0] + "-" + ref_name[1][0].split("_")[0]
            
        else:
            self.start_kmers = set([ref_seq[0:jf.k]])
            self.end_kmers = set([ref_seq[-jf.k:]])
            #if type(ref_seq) == str:
            #    ref_seq = [ref_seq]
            #    ref_name = [ref_name]
            #self.start_kmers = set([seq[0:jf.k] for seq in ref_seq])
            #self.end_kmers = set([seq[-jf.k:] for seq in ref_seq])
            
            #self.ref_mers = [[self.first_kmer] + uc.get_ref_kmers(seq, jf.k, name) + [self.last_kmer]
            #                 for seq, name in zip(ref_seq, ref_name)]
            #self.ref_set = set([kmer for seq in self.ref_mers for kmer in seq])
            self.ref_mers = [[self.first_kmer] + uc.get_ref_kmers(ref_seq, jf.k, ref_name) +\
                             [self.last_kmer]]
            self.ref_set = set(self.ref_mers[0])
        
        log.debug("Ref. set contains %d kmers.", len(self.ref_set))
        
        self.ref_seq = ref_seq
        self.ref_name = ref_name
        self.jf = jf
        self.mode = mode
        self.node_data = {}
        self.done = set()
        
        self.max_stack = max_stack
        self.max_break = max_break
        
        # register all k-mers from the ref
        for s in self.ref_set:
            if not (s == self.first_kmer or s == self.last_kmer):
                self.node_data[s] = self.jf.query(s)
        
        # kmer walking from each k-mer of ref_seq
        self.done.update(self.ref_set)
        for kmer in self.ref_set:
            if kmer in self.end_kmers or kmer == self.first_kmer or kmer == self.last_kmer:
                continue
            self.__extend([kmer], 0)
        
        self.graph_analysis(graphical)
    
    
    def __extend(self, stack, breaks):
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
                for p in stack:
                    self.node_data[p] = self.jf.query(p)
            else:
                self.__extend(stack + [child], breaks)
    
    
    def graph_analysis(self, graphical=False):
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
        
        self.paths = [] # in case there aren't any
        kmer = self.node_data.keys()
        kmer.extend([self.first_kmer, self.last_kmer])
        
        self.first_kmer_index = kmer.index(self.first_kmer)  # These will be useful after they are
        self.last_kmer_index = kmer.index(self.last_kmer)    # deleted from the kmer list
        
        num_k = len(kmer)
        graph = ug.Graph(num_k)
        
        log.debug("k-mer graph contains %d nodes.", num_k)
                
        log.debug("BigBang=%d, BigCrunch=%d" % (self.first_kmer_index, self.last_kmer_index))
        for s in self.start_kmers:
            log.debug("Start kmer %d %s" % (kmer.index(s), s))
        for e in self.end_kmers:
            log.debug("End   kmer %d %s" % (kmer.index(e), e))
         
        # The reference path, with node numbers
        ref_index = []
        # Look up indexes of chopped up k-mers from individual sequences in the nodes graph
        for i in range(len(self.ref_mers)):
            ref_index.append(map(lambda k: kmer.index(k), self.ref_mers[i]))
        
        # Build graph by finding pairwise kmer continuation from nodes
        for i in range(num_k):
            for j in range(num_k):
                if i == j:
                    continue
                if kmer[i][1:] == kmer[j][:-1]:
                    weight = 1
                    graph[i, j] = weight
        
        # Attribute a weight of 0.001 to continuous k-mers in the same sequence
        for l in range(len(ref_index)):  # for each sequence
            for k in range(len(ref_index[l])-1):  # k-mers in sequence - 1
                i = ref_index[l][k]
                j = ref_index[l][k+1]
                if i == self.first_kmer_index:
                    if sum([j in r for r in ref_index]) > 1:
                        continue
                graph[i, j] = 0.001
                # NOTE: A weight difference fold of 1000x might start causing problems for
                #       deletions that are > 31,000 bp long
        
        graph.init_paths(kmer.index(self.first_kmer), kmer.index(self.last_kmer))
        short_paths = graph.all_shortest()
        #print "\n".join([str(a) for a in short_paths])
        short_paths = [p[1:-1] for p in short_paths]
        
        # Remove artificially added first and last kmers
        for i in range(len(self.ref_mers)):
            self.ref_mers[i] = [k for k in self.ref_mers[i]
                                if k != self.first_kmer and k != self.last_kmer]
        
        kmer = kmer[:-2]  # Get rid of those artificial kmers ASAP
        
        if self.mode == "fusion":
            # Recover reference sequences
            ref_index_seq = []
            ref_index_full = []
            fusion_names = []
            for sp in short_paths:
                left_side = [[y for y in x if y != self.first_kmer_index]
                             for x in ref_index if self.first_kmer_index in x]
                right_side = [[y for y in x if y != self.last_kmer_index]
                              for x in ref_index if self.last_kmer_index in x]
                
                exon_num_left = set(range(len(left_side)))
                exon_num_right = set(range(len(right_side)))
                for i, k in enumerate(sp):
                    if exon_num_left:
                        bad_references = set()
                        for exon_num in exon_num_left:
                            target = left_side[exon_num]
                            if i < len(target) and k == target[i]:
                                continue
                            else:
                                bad_references.add(exon_num)
                        if len(exon_num_left) - len(bad_references):
                            exon_num_left = exon_num_left - bad_references
                            continue
                        else:
                            exon_num_left = set([sorted(list(exon_num_left),
                                                key=lambda x: len(left_side[x]))[0]])
                            exon_num = exon_num_left.pop()
                            ref_exon = get_seq(left_side[exon_num], kmer, skip_prefix=False)
                            fusion_names.append(self.ref_name_fusion[0][exon_num])
                            ref_index_seq.append(ref_exon)
                
                for i, k in enumerate(sp[::-1]):
                    if exon_num_right:
                        bad_references = set()
                        for exon_num in exon_num_right:
                            target = right_side[exon_num][::-1]
                            if i < len(target) and k == target[i]:
                                continue
                            else:
                                bad_references.add(exon_num)
                        if len(exon_num_right) - len(bad_references):
                            exon_num_right = exon_num_right - bad_references
                            continue
                        else:
                            exon_num_right = set([sorted(list(exon_num_right),
                                                        key=lambda x: len(right_side[x]))[0]])
                            exon_num = exon_num_right.pop()
                            ref_exon = get_seq(right_side[exon_num], kmer, skip_prefix=False)
                            fusion_names[-1] += "::" + self.ref_name_fusion[1][exon_num]
                            ref_index_seq[-1] += ref_exon
            
            ref_index_kmer = []
            for seq in ref_index_seq:
                ref_fusion_kmer = uc.get_ref_kmers(seq, self.jf.k, "ref_fusion")
                ref_index_kmer.append(ref_fusion_kmer)
                for s in ref_fusion_kmer:
                    if s not in kmer:  # slow (iterating through a list)
                        self.node_data[s] = self.jf.query(s)
                        kmer.append(s)
            
            for seq in ref_index_kmer:
                ref_index_full.append(map(lambda k: kmer.index(k), seq))
            
            ref_index = []  # just in case
            
        else:
            for i in range(len(ref_index)):
                ref_index[i] = [k for k in ref_index[i]
                                if k != self.first_kmer_index and k != self.last_kmer_index]
            ref_index = ref_index[0]  # temporary fix
        
        
        def get_name(a, b, offset=0):
            if self.mode == "fusion":
                try:
                    fusion = fusion_names[ref_index_full.index(a)]
                    return "Fusion\t{}".format(fusion)
                except ValueError:
                    ""  # We're in cluster mode
            
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
            for ind, path in enumerate(short_paths):
                if self.mode == "fusion":
                    ref_index = ref_index_full[ind]
                quant = upq.PathQuant(all_path=[path, ref_index],
                                      counts=self.node_data.values())
                
                quant.compute_coef()
                quant.refine_coef()
                quant.get_ratio()
                
                # Reference
                if list(path) == ref_index and self.mode != "fusion":
                    quant.adjust_for_reference()
                
                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name(ref_index, path),
                    seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                    ref_path=ref_index, info="vs_ref",
                    get_min_f=lambda path: min(get_counts(path, kmer)))
                
                self.paths += self.paths_quant
            
            if graphical:
                import matplotlib.pyplot as plt
                
                plt.figure(figsize=(10, 6))
                for path in short_paths:
                    plt.plot(get_counts(path, kmer),
                             label=get_name(ref_index, path).replace("\t", " "))
                plt.legend()
                plt.show()
        
        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        cluster = True
        if cluster:
            variant_diffs = []
            variant_set = set(range(0, len(short_paths)))
            for ind, variant in enumerate(short_paths):
                if self.mode == "fusion":
                    ref_index = ref_index_full[ind]
                diff = graph.diff_path_without_overlap(
                    ref_index, variant, self.jf.k)
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
                if self.mode == "fusion":
                    ref_index = ref_index_full[var_gr[2][0]]
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

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name(ref_path, path, offset),
                    seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                    ref_path=ref_path,
                    info="cluster %d n=%d" % (num_cluster, len(var_gr[2])),
                    get_min_f=lambda path: min(get_counts(path, kmer)),
                    start_off=start_off)

                self.paths_quant

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

    def get_paths(self):
        return self.paths

    def get_paths_quant(self):
        return self.paths_quant

    @staticmethod
    def output_header():
        upq.PathQuant.output_header()
