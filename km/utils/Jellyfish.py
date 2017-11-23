#                             -*- Mode: Python -*-
# Jellyfish.py --- Provides a python front-end to query a Jellyfish database.
#

import os
import sys
# import math
# import glob
# import re
# import subprocess
import logging as log
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

    def get_next(self, seq):
        # prev = self.query(seq)
        max_count = -1
        for c in ['A', 'C', 'G', 'T']:
            count = self.query(seq[1:] + c)
            # print c, ":", count
            if count > max_count:
                max_count = count
                next_c = c
        return (seq[1:] + next_c, max_count)

    def get_child(self, seq, forward=True):
        # count = self.query(seq)
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

        return map(lambda x: x[0], filter(lambda x: x[1] >= threshold, child))

    def get_best_child(self, seq, forward=True):
        count = self.query(seq)
        child = []
        sum = 0
        highest = 0
        for c in ['A', 'C', 'G', 'T']:
            if forward:
                c_seq = seq[1:] + c
            else:
                c_seq = c + seq[0:-1]
            c_count = self.query(c_seq)
            child.append((c_seq, c_count))
            if c_count > highest:
                highest = c_count
            sum += c_count
        threshold = max(sum * self.cutoff, self.n_cutoff)

        return map(lambda x: x[0], filter(lambda x: x[1] >= threshold and x[1] >= highest, child))

    def terminate(self):
        pass
        # self.jf.terminate()
