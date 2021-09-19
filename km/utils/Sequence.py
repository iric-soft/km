#                             -*- Mode: Python -*-
# Sequence.py
#

import sys

from . import common as uc


class RefSeq:
    """Reference path object that contains information on
    sequence, sequence name, kmers and kmer indices.

    Attributes
    ---------
    seq : str
        Reference sequence from fasta file.
    k : int
        K-mer size.
    first_kmer : str
        First kmer in seq.
    last_kmer : str
        Last kmer in seq.
    name : str
        Path name.
    ref_mer : list
        List of kmers found in seq.

    Methods
    -------
    set_index
    """

    def __init__(self, seq, name, k):

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


class AltSeq:
    """"""

    def __init__(self, alt_index, finder):
        """"""

        self.finder = finder

        self.seq_index = alt_index
        self.seq = self.finder.get_seq(self.seq_index, skip_prefix=False)
        self.first_ix = self.seq_index[0]
        self.last_ix = self.seq_index[-1]

        self.seq_len = len(self.seq_index)

        self.ref_index = self.finder.refpath.seq_index
        self.ref_name = self.finder.refpath.name

    def __getitem__(self, item):
        return self.seq_index[item]
