#                             -*- Mode: Python -*-
# Sequence.py
#

import sys

from . import common as uc


class RefSeq:
    """"""

    def __init__(self, seq, name, k):
        """"""""

        self.seq = seq
        self.k = k

        self.first_kmer = self.seq[0:self.k]
        self.last_kmer = self.seq[-self.k:]

        self.name = name

        self.ref_mer = uc.get_ref_kmer(seq, self.name, self.k)

        assert len(self.ref_mer)

    def set_index(self, kmer):
        """"""

        self.seq_index = tuple([kmer.index(k) for k in self.ref_mer])
        self.first_ix = self.seq_index[0]
        self.last_ix = self.seq_index[-1]

    def __getitem__(self, item):
        try:
            return self.seq_index[item]
        except AttributeError:
            sys.stderr.write('Attribute `seq_index` is not set yet\n')
            return None
