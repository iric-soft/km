#                             -*- Mode: Python -*-
# Jellyfish.py --- Provides a python front-end to query a Jellyfish database.
#

import os
import sys
import logging as log
try:
    import dna_jellyfish as jellyfish
except ModuleNotFoundError:
    import jellyfish

class Jellyfish:

    def __init__(self, filename, cutoff=0.30, n_cutoff=500, canonical=True):
        self.jf = jellyfish.QueryMerFile(filename)
        self.k = jellyfish.MerDNA.k()
        self.filename = filename
        self.cutoff = cutoff
        self.n_cutoff = n_cutoff
        self.canonical = canonical

    def query(self, seq):
        kmer = jellyfish.MerDNA(seq)
        if (self.canonical):
            kmer.canonicalize()
        return self.jf[kmer]

    def get_child(self, seq, forward=True):
        child = []
        sum = 0
        for c in ['A', 'C', 'G', 'T']:
            if forward:
                c_seq = seq[1:] + c
            else:
                c_seq = c + seq[0:-1]
            c_count = self.query(c_seq)
            child.append((c_seq, c_count))
            sum += c_count
        threshold = max(sum * self.cutoff, self.n_cutoff)

        return [x[0] for x in [x for x in child if x[1] >= threshold]]
