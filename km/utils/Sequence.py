#                             -*- Mode: Python -*-
# Sequence.py
#

import sys
import logging as log
from collections import namedtuple
from itertools import groupby

from . import common as uc


COUNTER = 0
NAMES = []

RefCandidate = namedtuple('RefCandidate',
    [
        'common',
        'divergence',
        'ref_index',
        'ref_name'
    ]
)


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

    def __init__(self, seq, attr, k):

        self.seq = seq
        self.attr = attr
        self.k = k

        self.first_kmer = self.seq[0:self.k]
        self.last_kmer = self.seq[-self.k:]

        self.__parse_attributes(attr)

        self.ref_mer = uc.get_ref_kmer(seq, self.name, self.k)
        assert len(self.ref_mer)

    def __parse_attributes(self, attr):
        try:
            exon_name = attr['gene_name']
        except KeyError:
            try:
                exon_name = attr['name']
            except KeyError:
                exon_name = attr['_filename']

        try:
            exon_id = attr["exon_id"]
        except KeyError:
            try:
                exon_id = attr["n"]
            except KeyError:
                #exon_id = hashlib.sha256(self.location.encode('utf-8'))
                #exon_id = exon_id.hexdigest()[:6].upper()
                global COUNTER
                COUNTER += 1
                exon_id = str(COUNTER)

        self.name = "%s.e%s" % (exon_name, exon_id)

        global NAMES
        assert self.name not in NAMES
        NAMES.append(self.name)

        exon_loc = attr["location"]
        if ":" not in exon_loc or "-" not in exon_loc:
            sys.exit(
                    'ERROR: Fasta entries do not contain a correctly ' +
                    'formatted location: {}\n'.format(exon_loc)
                )
        self.chromosome, loc = exon_loc.split(":")
        self.start, self.end = [int(x) for x in loc.split("-")]

        try:
            self.strand = attr["strand"]
        except KeyError:
            log.info("Attention! Strand is assumed to be '+' for exon %s.\n" % exon_id)
            self.strand = '+'

        self.location = '%s/%d-%d/%s' % (
                self.chromosome,
                self.start,
                self.end,
                self.strand
            )

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

    def __str__(self):
        return '%s=%s' % (
                self.name,
                self.location
            )


class BaseAltSeq:
    """"""

    def __init__(self, seq_index, finder):
        """"""

        self.finder = finder

        self.seq_index = seq_index
        self.seq = self.finder.get_seq(self.seq_index, skip_prefix=False)
        self.first_ix = self.seq_index[0]
        self.last_ix = self.seq_index[-1]

        self.seq_len = len(self.seq_index)

    def __getitem__(self, item):
        return self.seq_index[item]


class AltSeqSpawner(BaseAltSeq):
    def spawn(self):

        altpaths = list(self._spawn())

        # group paths by length
        paths_dct = {l: list(p) for l, p in groupby(altpaths, lambda x: x.seq_len)}
        lengths = sorted(list(paths_dct.keys()))[::-1]
        for pl in lengths:
            paths = paths_dct[pl]
            for i, path in enumerate(paths):
                skip = False

                # sanity check: paths cannot be the same length and equal
                for j in range(i+1, len(paths)):
                    path2 = paths[j]
                    assert path.seq != path2.seq

                longers = [x for x in lengths if x >= pl]
                while longers:
                    pl_longer = longers.pop()
                    for path2 in paths_dct[pl_longer]:
                        if path.seq in path2.seq and path.divergence > path2.divergence:
                            skip = True
                            longers = []
                            break
                if skip:
                    log.info('Skipping %s', path.ref_name)
                else:
                    yield path

    def _spawn(self):
        starts = [i for i, s in enumerate(self.seq_index) if s in self.finder.start_kmers_all_ix]
        ends = [j for j, e in enumerate(self.seq_index) if e in self.finder.end_kmers_all_ix]
        for i in starts:
            # include i positions for exons whose length is equal to self.jf.k
            new_ends = [j for j in ends if j > i]
            for j in new_ends:
                path = self.seq_index[i:j+1]
                nested_path = AltSeq(path, self.finder)
                nested_path.find_reference()
                yield nested_path


class AltSeq(BaseAltSeq):
    def find_reference(self):
        self.ref_index = None
        self.ref_name = None
        self.divergence = -1

        # Find all reference paths with start or end kmers that
        # correspond to this path
        start_candidates = [r for r in  self.finder.refpaths if r.first_ix == self.first_ix]
        end_candidates = [r for r in self.finder.refpaths if r.last_ix == self.last_ix]

        # Enumerate all potential pairs of references
        candidates = []
        for s in start_candidates:
            for e in end_candidates:
                if s.seq_index == e.seq_index:
                    candidates.append((s,))
                else:
                    candidates.append((s,e))

        candidates_mer = []
        for refs in candidates:
            if len(refs) == 1:
                ref_ix = refs[0].seq_index
                ref_name = refs[0].name
            else:
                junction = self.finder.kmer[refs[0][-1]] + self.finder.kmer[refs[1][0]]
                seq = uc.get_ref_kmer(junction, "junction", self.finder.jf.k)
                junction_mer = seq[1:-1]
                self.finder.add_kmers(junction_mer)
                junction_ix = tuple([self.finder.kmer.index(k) for k in junction_mer])

                ref_ix = refs[0].seq_index + junction_ix + refs[1].seq_index
                ref_name = refs[0].name + "::" + refs[1].name

            common = len(set(self.seq_index).intersection(set(ref_ix)))
            divergence = len(set(self.seq_index).union(set(ref_ix))) - common
            candidate = RefCandidate(common, divergence, ref_ix, ref_name)

            candidates_mer.append(candidate)

        highest_similarity = sorted(candidates_mer, key=lambda c: c.common, reverse=True)[0].common
        filtered_candidates = [c for c in candidates_mer if c.common == highest_similarity]
        if len(filtered_candidates) > 1:
            lowest_divergence = sorted(filtered_candidates, key=lambda c: c.divergence)[0].divergence
            filtered_candidates = [c for c in candidates_mer if c.divergence == lowest_divergence]
            if len(filtered_candidates) > 1:
                shortest_ref = len(sorted(filtered_candidates, key=lambda c: len(c.ref_index))[0].ref_index)
                filtered_candidates = [c for c in candidates_mer if len(c.ref_index) == shortest_ref]

        try:
            assert len(filtered_candidates) == 1
        except AssertionError:
            for fc in filtered_candidates:
                log.info(
                    "Warning: Candidate filtering failed. " +\
                    "Name=%s, Common=%d, Divergent=%d, RefLen=%d, RefSeq=%s, Seq=%s.",
                    fc.ref_name,
                    fc.common,
                    fc.divergence,
                    len(fc.ref_index),
                    self.finder.get_seq(fc.ref_index, False),
                    self.finder.get_seq(self.seq_index, False)
                )
        winner = filtered_candidates[0]

        self.ref_index = winner.ref_index
        self.ref_name = winner.ref_name
        self.divergence = winner.divergence
