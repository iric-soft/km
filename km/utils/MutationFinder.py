#                             -*- Mode: Python -*-
# MutationFinder.py
#

import sys
import logging as log
from collections import namedtuple

from . import Graph as ug
from . import PathQuant as upq
from . import Sequence as us
from . import common as uc


sys.setrecursionlimit(10000)


PathDiff = namedtuple('PathDiff',
    [
        'start',
        'end_ref',
        'end_var',
        'kmers_ref',
        'kmers_var',
        'end_ref_overlap'
    ]
)


class MutationFinder:
    """This class is the core of km. It sets up the environment
    and intializes kmer data structures.

    The first step is performed automatically where the jellyfish
    database (JF DB) is queried for all kmers that fit the
    elongation criteria given as parameters.

    The following step generates the graph and deduces alternative
    paths or shortest paths.

    Finally, paths are quantified and stored in a format ready for
    printing. An optional step is to generate mutation clusters and
    quantify those too.

    Attributes
    ---------
    refpath : km.utils.Sequence.RefSeq
        Reference path object.
    first_kmer : str
        Source capping node.
    last_kmer : str
        Sink capping node.
    ref_set : set
        Set of all kmers found in ref_seq
    jf: km.utils.Jellyfish.Jellyfish
        Jellyfish object with an open file handle to the JF DB.
    max_stack : int
        Maximum kmer extensions before reference getback.
    max_break : int
        Maximum branches before reference getback.
    paths : list
        Alternative paths discovered by `PathQuant`.
    node_data : dict
        Maps kmer string sequences to their position in the JF DB.
    kmer : list
        List of all kmers fetched from JF DB.
    counts : list
        Count of each kmer queried from JF DB.
    num_k : int
        Count of all kmers fetched from JF DB.
    refpaths : list
        List of kmer indices for the reference sequence.

    Methods
    -------
    graph_analysis
    quantify_paths
    quantify_clusters
    diff_path_without_overlap
    get_seq
    get_name
    get_counts
    get_paths
    output_header
    """

    def __init__(self,
            refpath,
            jf,
            max_stack=500,
            max_break=10,
            max_node=10000
        ):

        self.refpath = refpath

        self.first_seq = "BigBang"
        self.last_seq = "BigCrunch"

        self.ref_set = set(self.refpath.ref_mer)
        log.info("Ref. set contains %d kmers.", len(self.ref_set))

        self.jf = jf
        self.max_stack = max_stack
        self.max_break = max_break
        self.max_node = max_node

        self.node_data = {}

        # register all k-mers from the ref
        for s in self.ref_set:
            self.node_data[s] = self.jf.query(s)

        # kmer walking from each k-mer of ref_seq
        for seq in self.ref_set:
            # note: rightmost exons will extend purposelessly
            # because we removed check:
            # if seq == self.last_seq:
            #     continue
            self.__extend([seq])

        self.kmer = list(self.node_data.keys()) + [self.first_seq, self.last_seq]
        self.counts = list(self.node_data.values()) + [-1, -1]
        self.num_k = len(self.kmer)

        log.info("k-mer graph contains %d nodes.", self.num_k)

        # reference path, with node indices
        self.refpath.set_index(self.kmer)

        self.first_seq_ix = self.kmer.index(self.first_seq)
        self.last_seq_ix = self.kmer.index(self.last_seq)

        self.__define_edges()


    def __extend(self, stack, breaks=0):
        """Recursive depth first search"""

        if len(stack) > self.max_stack:
            return

        if len(self.node_data.keys()) > self.max_node:
            sys.exit(
                "ERROR: Node query count limit exceeded: max={}".format(
                    self.max_node
                )
            )

        cur_seq = stack[-1]
        childs = self.jf.get_child(cur_seq, forward=True)

        if len(childs) > 1:
            breaks += 1
            if breaks > self.max_break:
                return

        for child in childs:
            if child in self.node_data or child in set(stack):
                if child in set(stack) and child not in self.node_data:
                    log.info('Broke loop at kmer: %s' % child)
                for p in stack:
                    self.node_data[p] = self.jf.query(p)
            else:
                self.__extend(stack + [child], breaks)


    def __define_edges(self):
        self.start_kmers = set()
        self.end_kmers = set()

        if type(self.refpath) == us.RefSeq:
            first = self.refpath.first_kmer
            last = self.refpath.last_kmer

            self.start_kmers.add(first)

            self.end_kmers.add(last)

        self.start_kmers_ix = set([self.kmer.index(k) for k in self.start_kmers])
        self.end_kmers_ix = set([self.kmer.index(k) for k in self.end_kmers])

        log.info("BigBang=%d, BigCrunch=%d" % (self.first_seq_ix, self.last_seq_ix))
        for s in self.start_kmers:
            log.info("Start kmer %s %d" % (s, self.kmer.index(s)))
        for e in self.end_kmers:
            log.info("End   kmer %s %d" % (e, self.kmer.index(e)))


    @staticmethod
    def diff_path_without_overlap(ref, seq, k):
        """Compare a reference sequence to an alternative path
        discovered through kmer walking.

        Results from this function are used to pinpoint the exact
        location of the occurance of a mutation as well as the type
        of mutation detected.

        Here are some illustrative examples for the most common use-cases.

        - Reference (k=3)
            ref
                | ``(•••)•••••••••••••(•••)``
                | ``(000)0000000111111(1  )``
                | ``(0  )3456789012345(6  )``
            alt
                | ``(•••)•••••••••••••(•••)``
                | ``(000)0000000111111(1  )``
                | ``(0  )3456789012345(6  )``
            calculcated positions
                | ``> i     = 16``
                | ``> j_ref = 16``
                | ``> j_alt = 16``
                | ``> k_ref = 16 (unused)``

        - Substitution (k=3)
            ref
                | ``•••••••(•••)(•••)•••(•••)``
                | ``0000000(000)(111)111(1  )``
                | ``0123456(7  )(0  )345(6  )``
            alt
                | ``•••••••(••*)(•••)•••(•••)``
                | ``0000000(000)(111)111(1  )``
                | ``0123456(7  )(0  )345(6  )``
            calculcated positions
                | ``> i     = 07``
                | ``> j_ref = 10``
                | ``> j_alt = 10``
                | ``> k_ref = 10 (unused)``

        - ITD (k=3)
            ref
                | ``(•••)•••••••••••••[•••]••(•••)``
                | ``(000)0000000111111[1  ]12(2  )``
                | ``(0  )3456789012345[6  ]90(1  )``
            alt
                | ``(•••)•••••••••••••[•••]•••••••••••••[•••]••(•••)``
                | ``(000)0000000111111[111]1222222222233[333]33(3  )``
                | ``(0  )3456789012345[6  ]9012345678901[2  ]56(7  )``
            calculcated positions
                | ``> i     = 16``
                | ``> j_ref = 19``
                | ``> j_alt = 35``
                | ``> k_ref = 16``

        - Terminal ITD (k=3)
            ref
                | ``(•••)•••••••••••••[•••]``
                | ``(000)0000000111111[1  ]``
                | ``(0  )3456789012345[6  ]``
            alt
                | ``(•••)•••••••••••••[•••]•••••••••••••[•••]``
                | ``(000)0000000111111[111]1222222222233[3  ]``
                | ``(0  )3456789012345[6  ]9012345678901[2  ]``
            calculcated positions
                | ``> i     = 16``
                | ``> j_ref = 16``
                | ``> j_alt = 32``
                | ``> k_ref = 16``

        - Insertion (k=3)
            ref
                | ``• • • • • •(• • •)       (• • •)• • • •(• • •)``
                | ``0 0 0 0 0 0(0 0 0)       (0 1 1)1 1 1 1(1    )``
                | ``0 1 2 3 4 5(6    )       (9    )2 3 4 5(6    )``
            alt
                | ``• • • • • • •(• • *)* * *(• • •)• • • •(• • •)``
                | ``0 0 0 0 0 0 0(0 0 0)1 1 1(1 1 1)1 1 1 1(2    )``
                | ``0 1 2 3 4 5 6(7    )0 1 2(3    )6 7 8 9(0    )``
            calculcated positions
                | ``> i     = 07``
                | ``> j_ref = 09``
                | ``> j_alt = 13``
                | ``> k_ref = 09 (unused)``

        - Deletion (k=3)
            ref
                | ``• • • • • • •(• • •)• • •(• • •)•(• • •)``
                | ``0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1(1    )``
                | ``0 1 2 3 4 5 6(7    )0 1 2(3    )6(7    )``
            alt
                | ``• • • • • •(• • •)       (• • •)•(• • •)``
                | ``0 0 0 0 0 0 0 0 0         0 1 1 1(1    )``
                | ``0 1 2 3 4 5(6    )       (9    )2(3    )``
            calculcated positions
                | ``> i     = 07``
                | ``> j_ref = 13``
                | ``> j_alt = 09``
                | ``> k_ref = 13 (unused)``

        Parameters
        ---------
        ref : list
            List of kmer indices for the reference sequence.
        seq : list
            List of kmer indices for the alternative sequence.

        Returns
        -------
        namedtuple
            Named tuple with the following fields
                - start
                    Mutation start position (inclusive)
                - end_ref
                    Mutation end position on reference (last
                    position or last mutated kmer + 1)
                - end_var
                    Mutation end position on variant sequence (last
                    position or last mutated kmer + 1)
                - kmers_ref
                    List of kmer indices specific to the reference /
                    Deletion
                - kmers_var
                    List of kmer indices specific to the variant
                    sequence / Insertion
                - end_ref_overlap
                    Mutation end position on reference without k-mer
                    overlap; only used to detect ITDs for now
        """

        i = 0
        # Look for first left-side mutated kmer
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
        # Look for:
        #   - first right-side mutated kmer (if indel)
        # or
        #   - i + 31 (if snp or mnp)
        #
        # Conditions to look out for:
        #   - if insertion: j_seq > i + (k - 1)
        #   - if deletion: j_ref > i + (k - 1)
        #
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
        # Look for first right-side mutated kmer in ref (ignoring 'k' -> with overlap)
        # Useful for ITD detection
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
        """Generate a string from a list of kmer indices.

        Parameters
        ----------
        path : list
            The list of kmer IDs to join.
        skip_prefix : bool, optional
            If set to True only keep last nt. from first kmer.

        Returns
        -------
        str
            The resulting merge of kmers given as input.
        """

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
        """Compare a reference path to an alternative path and deduce
        mutation type and positions.

        Parameters
        ----------
        ref_ix : list
            Kmer IDs for the reference.
        path_ix : list
            Kmer IDs for the alternative path.
        offset : int, optional
            A precalculated offset measure to include in the output.

        Raises
        ------
        ...

        Returns
        -------
        str
            A string containing the mutation type and its positions
            in the reference sequence.
        """

        diff = self.diff_path_without_overlap(ref_ix, path_ix, self.jf.k)

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

        # Trim common end sequence (right-side) between del and ins
        trim = 1  # cannot be 0 because we use inverse indexing
        if len(del_seq) > 0:
            assert del_seq != ins_seq  # should never happen
            while del_seq[-trim:] == ins_seq[-trim:]:
                trim += 1
        trim -= 1  # offset by 1 to use as a right-side indexing extremity

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
        """Return kmer counts fetched from JF DB."""

        counts = [self.node_data[self.kmer[k]] for k in path]
        return counts

    def graph_analysis(self):
        """Perform kmer walking and find alternative paths

        Generate a 2-dimensional graph with all kmers queried
        at initialization. Nodes represent individual kmers and
        weighted edges are used to do the walking by prioritizing
        alternative paths.
        """

        self.paths = []

        # Initialize graph
        graph = ug.Graph(self.num_k)

        # For all kmers, find the next kmers with k - 1 overlap
        # and assign a weight of 1
        weight = 1

        # Match contiguous kmers
        prefix_dct = {}
        for i in range(self.num_k):
            prefix = self.kmer[i][:-1]
            try:
                prefix_dct[prefix].add(i)
            except KeyError:
                prefix_dct[prefix] = set([i])

        for i in range(self.num_k):
            suffix = self.kmer[i][1:]
            try:
                matches = prefix_dct[suffix]
            except KeyError:
                matches = set()
            for j in matches:
                if i != j:  # useful?
                    graph[i, j] = weight

        # Change the weights for contiguous kmers in reference
        # from 1 to 0.01 in graph
        weight = 0.01

        def adjust_graph_weights(ref_index):
            for k in range(len(ref_index)-1):
                i = ref_index[k]
                j = ref_index[k+1]
                graph[i, j] = weight

        adjust_graph_weights(self.refpath.seq_index)

        first_ix = self.kmer.index(self.first_seq)
        for start_ix in self.start_kmers_ix:
            graph[first_ix, start_ix] = weight

        last_ix = self.kmer.index(self.last_seq)
        for end_ix in self.end_kmers_ix:
            graph[end_ix, last_ix] = weight

        # Initialize paths and remove reference edges
        graph.init_paths(
            self.kmer.index(self.first_seq),
            self.kmer.index(self.last_seq)
        )

        # Locate shortest paths from non-reference edges
        short_paths = graph.all_shortest()
        # remove capping nodes from paths
        short_paths = [tuple(p[1:-1]) for p in short_paths]
        self.alt_paths = [us.AltSeq(s, self) for s in short_paths]

        # Group alternative paths with the same origin (ref_name) together
        alt_groups = {}
        for path in self.alt_paths:
            try:
                alt_groups[path.ref_name].append(path)
            except KeyError:
                alt_groups[path.ref_name] = [path]
        self.alt_groups = alt_groups


    def quantify_paths(self, graphical=False):
        """Quantify paths independently.

        Go through shortest paths one by one and quantify their
        expression. The result from the quantification method
        will be used in the final output.

        After quantification, paths are formatted and compiled
        into their final format for printing.

        Parameters
        ----------
        graphical : bool, optional
            If True generate a plot showing kmer coverage.
        """

        if graphical:
            import matplotlib.pyplot as plt

        def plot(paths):
            plt.figure(figsize=(10, 6))
            for path in paths:
                ref_name, ref_index, alt_index = path.ref_name, path.ref_index, path.seq_index
                plt.plot(
                    self.get_counts(alt_index),
                    label=self.get_name(
                            ref_index,
                            alt_index
                        ).replace("\t", " ") +\
                        ' (%s)' % ref_name
                )
            plt.legend()
            plt.show()

        if graphical:
            for paths in self.alt_groups.values():
                plot(paths)

        for path in self.alt_paths:
            ref_name, ref_index, alt_index = path.ref_name, path.ref_index, path.seq_index

            log.info('Quantifying %s' % ref_name)

            quant = upq.PathQuant(
                all_paths=[alt_index, ref_index],
                counts=self.counts
            )
            quant.compute_coef()
            quant.refine_coef()
            quant.get_ratio()

            # Reference
            if alt_index == ref_index:
                quant.adjust_for_reference()

            rvaf, ref_rvaf = quant.rVAF
            coef, ref_coef = quant.coef

            path_o = upq.Path(
                self.jf.filename,
                ref_name,
                self.get_name(ref_index, alt_index),
                rvaf,
                coef,
                min(self.get_counts(alt_index)),
                0,
                self.get_seq(alt_index, skip_prefix=False),
                ref_rvaf,
                ref_coef,
                self.get_seq(ref_index, skip_prefix=False),
                "vs_ref"
            )

            self.paths.append(path_o)


    def _find_clusters(self, alt_paths):
        """Generate clusters by cutting the sequence around mutations
        considering overlapping mutations as a cluster.
        """

        variant_diffs = []
        variant_set = set(range(0, len(alt_paths)))
        for path in alt_paths:
            ref_name, ref_index, alt_index = path.ref_name, path.ref_index, path.seq_index
            diff = self.diff_path_without_overlap(
                ref_index, alt_index, self.jf.k
            )
            variant_diffs.append(diff)

        def get_intersect(start, stop):
            for var in variant_set:
                cur_start = variant_diffs[var].start
                cur_end = variant_diffs[var].end_ref
                if cur_end >= start and cur_start <= stop:
                    if start == stop == cur_start == cur_end:
                        log.info('Terminal ITD ignored in cluster mode.')
                    elif stop == cur_end and (start == stop or cur_start == cur_end):
                        # meaning one is the reference and the other ends at ref end,
                        # which can happen when the ITD instantly checks at end
                        # extremities because it ends at `< end - k`
                        log.info('Quasi-terminal ITD ignored in cluster mode.')
                    else:
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

        # we know they all share the same reference
        ref_index = alt_paths[0].ref_index
        ref_name = alt_paths[0].ref_name

        for var_gr in variant_groups:
            start, stop, grp_ixs = var_gr

            if len(grp_ixs) == 1:
                var = grp_ixs[0]
                path = alt_paths[var]
                if path.seq_index == path.ref_index:
                    continue

            var_diffs = [variant_diffs[v] for v in grp_ixs]
            var_size = max(
                [abs(d.end_var - d.end_ref + 1) for d in var_diffs]
            )
            offset = max(0, start - var_size)
            ref_path = tuple(ref_index[offset:stop])

            clipped_paths = []
            for var in grp_ixs:
                cur_diff = variant_diffs[var]
                stop_off = cur_diff.end_var + stop - cur_diff.end_ref
                new_path = tuple(alt_paths[var][offset:stop_off])
                clipped_paths.append(new_path)

            yield (ref_name, ref_path, clipped_paths, offset)


    def quantify_clusters(self, graphical=False):
        """Detect and quantify cluster groups.

        In some cases, the complete sequence will contain at least
        2 homozygous mutations, which causes the overall minimumn
        coverage to be 0 for all paths. By defining clusters, we
        only keep the part of the sequence where there is only one
        mutation. Additinally, when multiple mutations overlap, they
        are grouped together into a cluster in order to get
        quantified as a single group. Then, go through clusters and
        quantify them in groups. The rest is similar to
        `quantify_paths`.

        After quantification, paths are formatted and compiled
        into their final format for printing.

        Parameters
        ----------
        graphical : bool, optional
            If True generate a plot showing kmer coverage for
            each cluster.
        """

        clusters = []
        for common_ref in sorted(self.alt_groups.keys(), key=lambda x: uc.natsortkey(x)):
            alt_paths = self.alt_groups[common_ref]
            for cluster in self._find_clusters(alt_paths):
                clusters.append(cluster)

        if graphical:
            import matplotlib.pyplot as plt

        for i, cluster in enumerate(clusters):
            ref_name, ref_path, clipped_paths, start_off = cluster
            num_cluster = i + 1

            log.info('Quantifying %s in cluster mode' % ref_name)

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
                counts=self.counts
            )
            quant.compute_coef()
            quant.refine_coef()
            quant.get_ratio()

            ref_rvaf, paths_rvaf = quant.rVAF[0], quant.rVAF[1:]
            ref_coef, paths_coef = quant.coef[0], quant.coef[1:]

            for path, rvaf, coef in zip(clipped_paths, paths_rvaf, paths_coef):
                assert path != ref_path
                path_o = upq.Path(
                    self.jf.filename,
                    ref_name,
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
        """Return all quantified paths from kmer walking.

        To be run after path quantification.

        Parameters
        ----------
        sort : bool, optional
            If True, sort based on variant position, type and
            minimum coverage.
        """

        if sort:
            paths = sorted(
                self.paths,
                key=lambda x: uc.natsortkey(*x[11].split(' '), x[1], x[3], x[2], x[6], rev_ix=[0])
            )
        else:
            paths = self.paths

        return paths

    @staticmethod
    def output_header():
        """Prints header output."""

        upq.Path.output_header()
