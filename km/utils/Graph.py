#                             -*- Mode: Python -*-
# Graph.py --- Generic graph class (directed, weighted)
#

import logging as log
import numpy as np

class Graph:
    """Nodes are linked with weighted edges dictating the shortest
    paths from each edge towards start and end nodes.

    Nodes are identified by their index: ``0 .. (n-1)``.

    Attributes
    ---------
    n : int
        Number of nodes (N).
    w : (N,N) array-like
        Graph with N nodes.
    edge_set : set
        A set of all relevant edges.
    first_node : int
        First node index.
    last_node : int
        Last node index.
    before : (N,) array-like
        The value at a specific position inside the list represents
        the node ID *preceding* the node whose ID is the current
        position.
    after : (N,) array-like
        The value at a specific position inside the list represents
        the node ID *following* the node whose ID is the current
        position.

    Methods
    -------
    init_paths
    all_shortest
    """

    def __init__(self, n):
        self.n = n
        self.w = np.empty((n, n), dtype=np.float32)
        self.w.fill(np.inf)
        self.edge_set = set()

        self.first_node = 0
        self.last_node = 0
        self.before = None
        self.after = None

        #self.ref_path = []

    def __getitem__(self, indices):
        (i, j) = indices
        return self.w[i, j]

    def __setitem__(self, indices, new_val):
        (i, j) = indices
        self.w[i, j] = new_val
        self.edge_set.add((i, j))

    def _get_paths(self, start, w):
        """Follow path from start until exhaustion of all kmers.

        Arguments
        ---------
        start : int
            ID of initial node to start with.
        w : numpy.ndarray
            A 2-dimensional weighted array containing weight
            values for each pair of edges.

        Returns
        -------
        numpy.ndarray
            List of preceding nodes as dictated by the weights
        """

        # a list of kmers (indices) and their corresponding
        # predecessors (values)
        prev = np.empty(self.n, dtype=np.int32)
        prev.fill(-1)

        # dist will contain the cumulative weights of forward
        # kmers as they appear during kmer walking
        dist = np.empty(self.n, dtype=np.float32)
        dist.fill(np.inf)
        unvisited = np.ones(self.n, dtype=bool)

        def visit(i):
            # ndist is the weights for all kmers forward of `i` +
            # the current cumulative weight at the stage of kmer `i`
            ndist = w[i, :] + dist[i]

            # only keep forward kmers which:
            #   - were not selected before (np.inf in dist)
            #   - are currently in a reference edge even if they
            #     were already found in a variant edge, their
            #     value in dist will be overwritten by the lower
            #     weight and their previous kmer will be the current
            #     kmer `i`.
            ind = tuple([ndist < dist])
            dist[ind] = ndist[ind]
            prev[ind] = i

        # ensure we start with the start kmer, at this point dist
        # is filled with +np.inf except for `start` which is 0
        dist[start] = 0
        while any(unvisited):
            # select passed over kmer with smallest weight from
            # previous iterations
            unv_ix = np.where(unvisited)[0]
            i = unv_ix[dist[unv_ix].argmin()]
            # find and annotate forward kmers for kmer `i`
            visit(i)
            unvisited[i] = False

        return prev

    def init_paths(self, first_node, last_node):
        """For each node, find its preceding (array `before`) and
        following (array `after`) node based on the weights in
        the graph.

        Construct a one-to-one mapping of nodes from each extremity
        of the reference, and get rid of reference edges (discovered
        by following the smallest weights from start to end).

        Follow the path of least-resistence / smallest weight
        starting from node `start` and continuously affect each
        node with its deduced predecessor.

        Follow paths from both sides of the reference to be able
        to construct the shortest paths from each side of an
        alternative edge.

        Here is a schematic representation of how the path is
        constructed. When 2 alternative branches exist, branched
        nodes are affected with the same preceding node. However,
        when the branching has to close, only one predecessor can
        be chosen, and it will always be the one with the lowest
        weight.

        | ``|------------------------------------------------------------------|``
        | ``|               • <- • <- •                                        |``
        | ``|             /             \   <= lower weight == is kept         |``
        | ``| • <- • <- •                 • <- • <- •                          |``
        | ``|             \            (/)  <= higher weight == is overwritten |``
        | ``|               • <- • <- •                                        |``
        | ``|------------------------------------------------------------------|``

        Notes
        -----
        In extremely rare situations, if the reference path (small
        weights) goes on for a very long time, the cumulative
        weights could end up catching up to the alternative path and
        become higher weights, which would cause an inversion
        of reference and alternative paths. An example would be if
        a reference path with a weight of 0.01 is 100 times longer
        than an alternative path with weights of 1, we may then get
        to a point where the cumulative weight of the reference
        is higher than that of the alternative.

        Arguments
        ---------
            first_node : int
                First node ID.
            last_node : int
                Last node ID.
        """

        self.first_node = first_node
        self.last_node = last_node
        self.before = self._get_paths(self.first_node, self.w)
        self.after = self._get_paths(self.last_node, self.w.transpose())

        # Load up and remove edges from the ref path. We know for a fact
        # that starting from the first node in list `after` will guarantee
        # that all edges will come from the reference (i.e., following the
        # path of least resistance). We would also get to the same result
        # by starting with the last node in list `before`.
        #path = [self.first_node]
        removed = 0
        curs = set(np.where(self.before == self.first_node)[0])
        for cur in curs:
            #path = [self.first_node]
            last_cur = None  # ensures we keep only one edge from the reference
            while self.after[cur] != -1:
                cur = self.after[cur]
                #path.append(cur)
                if last_cur and (last_cur, cur) in self.edge_set:  # we might have taken care of this edge already
                    self.edge_set.remove((last_cur, cur))
                    log.debug("Removing (%d, %d)", last_cur, cur)
                    removed += 1
                last_cur = cur
            #self.ref_path.append(path)
        log.info("Removed %d ref edges.", removed)

    def _get_shortest(self, a, b):
        """Return shortest path passing through edge (a, b)"""

        def follow(start, prev):
            cur = start
            while prev[cur] != -1:
                cur = prev[cur]
                path.append(cur)

        path = [b, a]
        follow(a, self.before)
        path.reverse()
        follow(b, self.after)

        ## Only keep paths from source to sink
        if path[0] != self.first_node or path[-1] != self.last_node:
            path = None

        return path

    def all_shortest(self):
        """Find every unique path that can be constructed from a
        non-reference edge by following the paths in `before` and
        `after`.

        Returns
        -------
        list
            List of unique shortest paths walked from all edges
        """

        log.info("%d edges in non-ref edge set.", len(self.edge_set))

        all_paths = set()
        for (i, j) in self.edge_set:
            log.debug("Computing shortest path through edge: (%d, %d)", i, j)
            path = self._get_shortest(i, j)
            if path:
                all_paths.add(tuple(path))

        return list(all_paths)
