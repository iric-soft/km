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
            self.ref_mers = [[self.first_kmer] +\
                             uc.get_ref_kmers(ref_seq, jf.k, ref_name) +\
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
                
        def get_name(a, b, offset=0):
            if self.mode == "fusion" and not cluster:
                fusion = fusion_names[ref_index_full.index(a)]
                return "Fusion\t{}".format(fusion)
            
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
        
        def locate_reference(path, sequences, exon_num, start_kmers):
            def find_all(path, exons_in_order, exon_alternatives):
                if not path:
                    exon_alternatives.append(exons_in_order)
                    return exon_alternatives
                exons_to_explore = []
                for i, k in enumerate(path):
                    if exons_to_explore:
                        #print exons_to_explore
                        for j, ee in enumerate(exons_to_explore):
                            if ee[4] or ee[1] == ee[2]:
                                ee[4] = True
                            elif k == sequences[ee[0]][ee[2] + ee[3]]:
                                ee[2] += 1
                                ee[2] += ee[3] # we had a snp or mnp
                                ee[3] = 0
                                ee[4] = False
                            elif k == sequences[ee[0]][-1]:
                                ee[4] = True
                            else:
                                ee[3] += 1
                                ee[4] = None
                        # If some are True and some are None get rid of the None
                        if sum([ee[4] for ee in exons_to_explore if ee[4]]):
                            #print "getting rid of None", exons_to_explore
                            exons_to_explore = [ee for ee in exons_to_explore if not ee[4] is None]
                        # If all are True
                        if sum([ee[4] for ee in exons_to_explore if ee[4]]) == len(exons_to_explore):
                            #print "all true", exons_to_explore
                            completed = [ee for ee in exons_to_explore if ee[4]]  # necessary?
                            #print "all true completed", completed
                            if exons_in_order:
                                exons_in_order.append(sorted(completed, key=lambda x: x[1])[-1][0])
                            else:
                                # Get all start alternatives from first exon and readjust paths later
                                if len(exons_to_explore) > 1:
                                    p = []
                                    for l in range(1, len(path)):
                                        if path[l] in start_kmers:
                                            p = path[l:]
                                            #print "ALTERNATIVE START"
                                            exon_alternatives.extend(locate_reference(
                                                p, sequences, exon_num, start_kmers))
                                            break
                                exons_in_order.append(sorted(completed, key=lambda x: x[1])[-1][0])
                            #print "exons in order", exons_in_order
                            p = []
                            for l in range(i, len(path)):
                                if path[l] in start_kmers:
                                    p = path[l:]
                                    break
                            return find_all(p, exons_in_order, exon_alternatives)
                    if k in start_kmers:
                        #print "START", k
                        for e in exon_num:
                            if sequences[e][0] == k:
                                exons_to_explore.append([e, len(sequences[e]), 1, 0, False])
            return find_all(path, [], [])
        
        cluster = False
        
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
                if i == self.first_kmer_index:  # remove extremities that are part of other paths
                    if sum([j in r[2:] for r in ref_index]):
                        continue
                elif j == self.last_kmer_index:
                    if sum([i in r[:-2] for r in ref_index]):
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
            # Recover reference sequences and correct paths to reference if needed
            ref_index_seq = []
            fusion_names = []
            new_paths = []
            for sp in short_paths:
                kleft = [[y for y in x if y != self.first_kmer_index]
                             for x in ref_index if self.first_kmer_index in x]
                kright = [[y for y in x if y != self.last_kmer_index]
                              for x in ref_index if self.last_kmer_index in x]
                kright_rev = [x[::-1] for x in kright]
                
                start_kmers_idx = map(lambda k: kmer.index(k), self.start_kmers)
                end_kmers_idx = map(lambda k: kmer.index(k), self.end_kmers)
                # indexes of exon candidates on left and right
                exon_num_left = range(len(kleft))
                exon_num_right = range(len(kright))
                
                #print "start"
                exons_start_all = locate_reference(sp, kleft, exon_num_left, start_kmers_idx)
                #print "end"
                exons_end_all = locate_reference(sp[::-1], kright_rev, exon_num_right, end_kmers_idx)
                
                #print "done"
                #print exons_start_all
                #print exons_end_all
                
                for exons_start in exons_start_all:
                    for exons_end in exons_end_all:
                        assert len(exons_end) == 1  # should never continue beyond the fusion exon
                        first_ex_k = kleft[exons_start[0]][0]
                        last_ex_k = kright[exons_end[0]][-1]
                        new_p = [s for s in sp]
                        for i, k in enumerate(new_p):
                            if k == first_ex_k:
                                new_p = new_p[i:]
                        for i, k in enumerate(new_p[::-1]):
                            if k == last_ex_k and i:
                                new_p = new_p[:-i]
                        if new_p in new_paths:
                            continue
                        new_paths.append(new_p)
                        ref_index_seq.append("")
                        fusion_names.append("")
                        for e in exons_start:
                            ref_exon = get_seq(kleft[e], kmer, skip_prefix=False)
                            ref_index_seq[-1] += ref_exon
                            fusion_names[-1] += self.ref_name_fusion[0][e] + "::"
                        ref_index_seq[-1] += get_seq(kright[exons_end[0]], kmer, skip_prefix=False)
                        fusion_names[-1] += self.ref_name_fusion[1][exons_end[0]]
                
            ref_index_kmer = []
            ref_index_full = []
            for seq in ref_index_seq:
                ref_fusion_kmer = uc.get_ref_kmers(seq, self.jf.k, "ref_fusion")
                ref_index_kmer.append(ref_fusion_kmer)
                for s in ref_fusion_kmer:
                    if s not in kmer:  # slow (iterating through a list)
                        self.node_data[s] = self.jf.query(s)
                        kmer.append(s)
            for seq in ref_index_kmer:
                ref_index_full.append(map(lambda k: kmer.index(k), seq))
            ref_index = ref_index_full
            short_paths = new_paths
        else:
            ref_index = [[x for x in ref_index[0] if x != self.first_kmer_index
                                                  and x != self.last_kmer_index]]*len(short_paths)
       
        # Quantify all paths independently
        individual = True
        if individual:
            for path, ref_ind in zip(short_paths, ref_index):
                quant = upq.PathQuant(all_path=[path, ref_ind],
                                      counts=self.node_data.values())
                
                quant.compute_coef()
                quant.refine_coef()
                quant.get_ratio()
                
                # Reference
                # TODO: Make it work with fusion. The way it is now returns 0.5 for the ratio
                if list(path) == ref_ind and self.mode != "fusion":
                    quant.adjust_for_reference()
                
                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name(ref_ind, path),
                    seq_f=lambda path: get_seq(path, kmer, skip_prefix=False),
                    ref_path=ref_ind, info="vs_ref",
                    get_min_f=lambda path: min(get_counts(path, kmer)))
                
                self.paths += self.paths_quant
            
            if graphical:
                import matplotlib.pyplot as plt
                
                plt.figure(figsize=(10, 6))
                for path in short_paths:
                    plt.plot(get_counts(path, kmer),
                             label=get_name(ref_ind, path).replace("\t", " "))
                plt.legend()
                plt.show()
        
        # Quantify by cutting the sequence around mutations,
        # considering overlapping mutations as a cluster
        cluster = True
        if cluster:
            variant_diffs = []
            variant_set = set(range(0, len(short_paths)))
            for variant, ref_ind in zip(short_paths, ref_index):
                diff = graph.diff_path_without_overlap(ref_ind, variant, self.jf.k)
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
                ref_ind = ref_index[var_gr[2][0]]
                if (len(var_gr[2]) == 1 and
                        list(short_paths[var_gr[2][0]]) == ref_ind):
                    continue
                num_cluster += 1

                start = var_gr[0]
                stop = var_gr[1]
                var_size = max([abs(x[2]-x[1]+1) for x in [variant_diffs[v] for v in var_gr[2]]])
                offset = max(0, start - var_size)
                ref_path = ref_ind[offset:stop]
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
