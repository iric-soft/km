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
                 max_break=10, max_node=10000, attr={}):
        
        # Load the reference sequence and prepare ref k-mers
        
        self.first_kmer = "BigBang"
        self.last_kmer = "BigCrunch"
        
        # Detect mode
        if attr and type(ref_seq) == list:
            genes = [at["name"] for at in attr]
            self.mode = "fusion"
        else:
            self.mode = "mutation"
        
        if self.mode == "mutation":
            self.start_kmers = set([ref_seq[0:jf.k]])
            self.end_kmers = set([ref_seq[-jf.k:]])
            
            self.ref_mers = [[self.first_kmer] +\
                             uc.get_ref_kmers(ref_seq, jf.k, ref_name) +\
                             [self.last_kmer]]
            
            self.ref_set = set(self.ref_mers[0])
            
        elif self.mode == "fusion":
            self.exons_name = ["%se%s" % (at["name"], at["n"]) for at in attr]
            
            self.start_kmers = set([seq[0:(jf.k)] for seq in ref_seq])
            self.end_kmers = set([seq[-(jf.k):] for seq in ref_seq])
            
            self.ref_mers = [[self.first_kmer] + uc.get_ref_kmers(seq, jf.k, name) + [self.last_kmer]
                             for seq, name in zip(ref_seq, self.exons_name)]
            self.ref_set = set([kmer for seq in self.ref_mers for kmer in seq])
        
        log.info("Ref. set contains %d kmers.", len(self.ref_set))
        
        self.ref_seq = ref_seq
        self.ref_name = ref_name
        self.jf = jf
        self.node_data = {}
        self.done = set()
        
        self.max_stack = max_stack
        self.max_break = max_break
        self.max_node = max_node
        
        # register all k-mers from the ref
        for s in self.ref_set:
            if not (s == self.first_kmer or s == self.last_kmer):
                self.node_data[s] = self.jf.query(s)
        
        # kmer walking from each k-mer of ref_seq
        self.done.update(self.ref_set)
        for kmer in self.ref_set:
            if kmer == self.first_kmer or kmer == self.last_kmer:
                # of note, rightmost exons will extend purposelessly
                continue
            self.__extend([kmer], 0)
        
        # TODO: Should eventually include that in all tools and remove from here
        if self.mode == "mutation":
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
                if len(self.done) > self.max_node:
                    sys.exit("ERROR: Node query count limit exceeded: max={}".format(self.max_node))
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
            if self.mode == "fusion":
                if cluster:
                    fn = fusion_names[ref_index.index(a[1])]
                    a = a[0]
                else:
                    fn = fusion_names[ref_index.index(a)]
            
                fusion = "/{}".format(fn)
            else:
                if cluster:
                    a = a[0]
                fusion = ""
            
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
                split_exons = fusion[1:].split("::")
                if len(split_exons) > 1:
                    if split_exons[0].rsplit("e")[0] != split_exons[1].rsplit("e")[0]:
                        return "Fusion{}\t".format(fusion)
                return "Reference{}\t".format(fusion)
            else:
                fus = ""
                split_exons = fusion[1:].split("::")
                if len(split_exons) > 1:
                    if split_exons[0].rsplit("e")[0] != split_exons[1].rsplit("e")[0]:
                        fus = "Fusion-"
                variant = "Indel"
                # SNP have equal length specific sequences
                if diff[1] == diff[2]:
                    variant = "Substitution"

                # ITD have zero kmers in ref after full trimming.
                # However, this does not distinguish cases where there is
                # garbage between repeats.
                elif diff[0] == diff[5]:
                    variant = "Insertion"  # "ITD"; Prone to errors, trust find_report instead
                elif len(del_seq) == 0 and len(ins_seq) != 0:
                    variant = "Insertion"
                elif len(del_seq) != 0 and len(ins_seq) == 0:
                    variant = "Deletion"

                if diff[0] + k > diff[1]:  # some special case of ITD
                    diff0 = diff[0] - (diff[0] + k - diff[1])
                else:
                    diff0 = diff[0]

                return "{}{}\t{}:{}:{}".format(
                    fus+variant,
                    fusion,
                    diff0 + k + offset,
                    (string.lower(del_seq) + "/" + ins_seq),
                    diff[1] + 1 + offset)
        
        def get_counts(path, kmer):
            counts = []
            for i in path:
                counts += [self.node_data[kmer[i]]]
            # print("length counts: " + str(len(counts)))
            # print("min counts: " + str(min(counts)))
            
            return counts
         
        cluster = False
        
        self.paths = [] # in case there aren't any
        kmer = self.node_data.keys()
        self.counts = self.node_data.values()
        kmer.extend([self.first_kmer, self.last_kmer])
        
        self.first_kmer_index = kmer.index(self.first_kmer)  # These will be useful after they are
        self.last_kmer_index = kmer.index(self.last_kmer)    # deleted from the kmer list
        
        num_k = len(kmer)
        graph = ug.Graph(num_k)
        
        log.info("k-mer graph contains %d nodes.", num_k)
         
        log.info("BigBang=%d, BigCrunch=%d" % (self.first_kmer_index, self.last_kmer_index))
        for s in self.start_kmers:
            log.info("Start kmer %d %s" % (kmer.index(s), s))
        for e in self.end_kmers:
            log.info("End   kmer %d %s" % (kmer.index(e), e))
        
        # The reference path, with node numbers
        ref_index = []
        # Look up indexes of chopped up k-mers from individual sequences in the nodes graph
        for i in range(len(self.ref_mers)):
            ref_index.append(map(lambda k: kmer.index(k), self.ref_mers[i]))
        
        # Build graph by finding pairwise kmer continuation from nodes
        kmer_start = [k[:-1] for k in kmer]
        kmer_end = [k[1:] for k in kmer]
        for i in range(num_k):
            starts = [j for j, k in enumerate(kmer_start) if k == kmer_end[i]]
            for j in starts:
                if i == j:
                    continue
                weight = 1
                graph[i, j] = weight
        #for i in range(num_k):
        #    for j in range(num_k):
        #        if i == j:
        #            continue
        #        if kmer[i][1:] == kmer[j][:-1]:
        #            weight = 1
        #            graph[i, j] = weight
        
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
        
        short_paths = [p[1:-1] for p in short_paths]
        
        # Remove artificially added first and last kmers
        for i in range(len(self.ref_mers)):
            self.ref_mers[i] = [k for k in self.ref_mers[i]
                                if k != self.first_kmer and k != self.last_kmer]
        kmer = kmer[:-2]  # Get rid of those artificial kmers ASAP
        ks = set(kmer)
        
        kall = [[y for y in x if y != self.first_kmer_index and y != self.last_kmer_index]
                     for x in ref_index]
        kall_index = {e:{k:i for i, k in enumerate(seq)} for e, seq in enumerate(kall)}
        
        if self.mode == "fusion":
            # Recover reference sequences and correct paths to reference if needed
            new_paths = []
            ref_index_seq = []
            ref_index_full = []
            fusion_names = []
            for sp in short_paths:
                exons_all = uc.locate_reference(
                        sp,
                        kall,
                        self.exons_name,
                        kall_index)
                 
                for start_exon, end_exon in exons_all: 
                    if end_exon is None:
                        start_exon_path = kall[start_exon]
                        first_ex_k = start_exon_path[0]
                        last_ex_k = start_exon_path[-1]
                        sp_ref_seq = get_seq(start_exon_path, kmer, False)
                        sp_ref_name = self.exons_name[start_exon]
                    else:
                        start_exon_path = kall[start_exon]
                        end_exon_path = kall[end_exon]
                        first_ex_k = start_exon_path[0]
                        last_ex_k = end_exon_path[-1]
                        sp_ref_seq = get_seq(start_exon_path, kmer, False) +\
                                     get_seq(end_exon_path, kmer, False)
                        sp_ref_name = self.exons_name[start_exon] + "::" + self.exons_name[end_exon]
                        # Append new reference indexes
                        first_junction_k = start_exon_path[-1:]
                        last_junction_k = end_exon_path[0:1]
                        junction_seq = get_seq(first_junction_k, kmer, False) +\
                                       get_seq(last_junction_k, kmer, False)
                        junction_kmer = uc.get_ref_kmers(junction_seq, self.jf.k, "junction_seq")
                        for k in junction_kmer:
                            if k not in ks:
                                ks.add(k)
                                c = self.jf.query(k)
                                kmer.append(k)
                                self.counts.append(c)
                                self.node_data[k] = c
                    
                    new_p = [s for s in sp]
                    for i, k in enumerate(sp):
                        if k == first_ex_k:
                            new_p = new_p[i:]
                            break  # for ITDs
                    for i, k in enumerate(sp[::-1]):
                        if k == last_ex_k and not i:  # equivalent to the non-working new_p[:-0]
                            break
                        if k == last_ex_k and i:
                            new_p = new_p[:-i]
                            break  # for ITDs
                    new_paths.append(new_p)
                    ref_index_seq.append(sp_ref_seq)
                    sp_ref_kmer = uc.get_ref_kmers(sp_ref_seq, self.jf.k, sp_ref_name, False)
                    ref_index_full.append(map(lambda k: kmer.index(k), sp_ref_kmer))
                    fusion_names.append(sp_ref_name)
            
            new_paths_filtered = []
            ref_index_full_filtered = []
            fusion_names_filtered = []
            for (i, p) in enumerate(new_paths):
                ref = ref_index_full[i]
                name = fusion_names[i]
                seq = get_seq(p, kmer, False)
                skip = False
                
                for (p2, ref2) in zip(new_paths_filtered, ref_index_full_filtered):
                    if p == p2 and ref == ref2:
                        log.info("Omitting duplicate path {}".format(name))
                        skip = True
                        break
                if skip:
                    continue
                
                for (j, p2) in enumerate(new_paths):
                    if i == j:
                        continue  # we're looking at the same path
                    ref2 = ref_index_full[j]
                    seq2 = get_seq(p2, kmer, False)
                    
                    if seq in seq2 and p != ref:  # not a reference
                        assert seq2 != ref2
                        diff1 = len(p) + len(ref) - len(set(p).intersection(set(ref)))*2
                        diff2 = len(p2) + len(ref2) - len(set(p2).intersection(set(ref2)))*2
                        if diff1 > diff2:  # Parcimony
                            if seq == seq2:
                                log.info("Omitting incorrect ref. {}".format(name))
                            else:
                                log.info("Omitting a nested path: {}".format(name))
                            skip = True
                            break
                if skip:
                    continue
                
                # If both loops say it's okay:
                new_paths_filtered.append(p)
                ref_index_full_filtered.append(ref)
                fusion_names_filtered.append(name)
             
            ref_index = ref_index_full_filtered
            short_paths = new_paths_filtered
            fusion_names = fusion_names_filtered
            #ref_index = ref_index_full
            #short_paths = new_paths
            #fusion_names = fusion_names
         
        elif self.mode == "mutation":
            ref_index = [[x for x in ref_index[0] if x != self.first_kmer_index
                                                  and x != self.last_kmer_index]]*len(short_paths)
        
        log.info("Number of paths found: {}".format(len(short_paths)))
       
        # Quantify all paths independently
        individual = True
        if individual:
            for path, ref_ind in zip(short_paths, ref_index):
                quant = upq.PathQuant(all_path=[path, ref_ind],
                                      counts=self.counts)
                
                quant.compute_coef()
                quant.refine_coef()
                quant.get_ratio()
                
                # Reference
                if list(path) == ref_ind:
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

            def get_intersect(start, stop, ref):
                for var in variant_set:
                    if (ref_index[var] == ref and
                            variant_diffs[var][1] >= start and
                            variant_diffs[var][0] <= stop):
                        return var
                return -1

            # CHANGES MADE FOR FUSION AND ISOFORMS MODE SHOULD BE CORRECTED
            variant_groups = []
            while len(variant_set) > 0:
                seed = variant_set.pop()
                ref_ind = ref_index[seed]
                grp = [seed]
                start = variant_diffs[seed][0]
                stop = variant_diffs[seed][1]
                variant = get_intersect(start, stop, ref_ind)
                while variant != -1:
                    variant_set.remove(variant)
                    grp += [variant]
                    start = min(start, variant_diffs[variant][0])
                    stop = max(stop, variant_diffs[variant][1])
                    variant = get_intersect(start, stop, ref_ind)
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
                                      counts=self.counts)

                quant.compute_coef()
                quant.refine_coef()

                quant.get_ratio()

                self.paths_quant = quant.get_paths(
                    db_f=self.jf.filename,
                    ref_name=self.ref_name,
                    name_f=lambda path: get_name([ref_path, ref_ind], path, offset),
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
                                     label=get_name([ref_path, ref_ind], path, offset).split("\t")[0])
                    plt.legend()
                    plt.show()

    def get_paths(self):
        return self.paths

    def get_paths_quant(self):
        return self.paths_quant

    @staticmethod
    def output_header():
        upq.PathQuant.output_header()
