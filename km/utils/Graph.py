#                             -*- Mode: Python -*-
# Graph.py --- Generic graph class (directed, weighted)
#

import logging as log
import numpy as np

class Graph:
    ## Nodes are identified by 0..(n-1)

    def __init__(self, n):
        self.n = n
        self.w = np.empty((n, n), dtype=np.float32)
        self.w.fill(np.inf)
        self.edge_set = set()

        self.first_node = 0
        self.last_node = 0
        self.before = None
        self.after = None

        self.ref_path = []

    def __getitem__(self, indices):
        (i, j) = indices
        return self.w[i, j]

    def __setitem__(self, indices, new_val):
        (i, j) = indices
        self.w[i, j] = new_val
        self.edge_set.add((i, j))

    def _get_paths(self, start, w):
        prev = np.empty(self.n, dtype=np.int32)
        prev.fill(-1)
        dist = np.empty(self.n, dtype=np.float32)
        dist.fill(np.inf)
        unvisited = set(range(self.n))

        def visit(i):
            ndist = w[i, :] + dist[i]
            for j in range(self.n):
                if ndist[j] < dist[j]:
                    dist[j] = ndist[j]
                    prev[j] = i

        dist[start] = 0
        visit(start)
        unvisited.remove(start)

        def min_unvisited():
            min_i = -1
            min_dist = np.inf
            for i in unvisited:
                if dist[i] <= min_dist:
                    min_i = i
                    min_dist = dist[i]
            return min_i

        while unvisited:
            i = min_unvisited()
            unvisited.remove(i)
            visit(i)

        return prev

    def init_paths(self, first_node, last_node):
        self.first_node = first_node
        self.last_node = last_node
        self.before = self._get_paths(self.first_node, self.w)
        self.after = self._get_paths(self.last_node, self.w.transpose())
        # Load up and remove edges from the ref path
        path = [self.first_node]
        cur = self.first_node
        last_cur = None
        while self.after[cur] != -1:
            cur = self.after[cur]
            path.append(cur)
            if last_cur:
                self.edge_set.remove((last_cur, cur))
                log.debug("Removing (%d, %d)", last_cur, cur)
            last_cur = cur
        self.ref_path = path

    def get_shortest(self, a, b):
        # Returns the shortest path passing through edge (a, b)
        path = [b, a]

        def follow(start, prev):
            cur = start
            while prev[cur] != -1:
                cur = prev[cur]
                path.append(cur)

        follow(a, self.before)
        path.reverse()
        follow(b, self.after)
        ## Only keep paths from source to sink
        if path[0] != self.first_node or path[-1] != self.last_node:
            return None
        return path

    def all_shortest(self):
        all_paths = set()
        log.debug("%d edges in non-ref edge set.", len(self.edge_set))
        for (i, j) in self.edge_set:
            log.debug("Computing shortest path through edge: (%d, %d)", i, j)
            path = self.get_shortest(i, j)
            if path:
                all_paths.add(tuple(path))
        return list(all_paths)

    def diff_path_without_overlap(self, ref, seq, k):
        # Returns (start, stop_ref, stop_variant, kmers_ref, kmers_variant, stop_ref_fully_trimmed)
        i = 0

        while i < len(ref) and i < len(seq) and ref[i] == seq[i]:
            i += 1

        j_ref = len(ref)
        j_seq = len(seq)
        while j_ref > i + (k - 1) and j_seq > i + (k - 1) and ref[j_ref - 1] == seq[j_seq - 1]: #  + (k - 1):to prevent kmer from overlapping
            j_ref -= 1
            j_seq -= 1

        k_ref = j_ref
        k_seq = j_seq
        while k_ref > i and ref[k_ref - 1] == seq[k_seq - 1]:
            k_ref -= 1
            k_seq -= 1

        # log.debug("diffpath : " + " ".join (str(x) for x in [i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref]))

        return (i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref)
