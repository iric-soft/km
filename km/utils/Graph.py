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

        #self.ref_path = []

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
        unvisited = np.ones(self.n, dtype=bool)

        def visit(i):
            ndist = w[i, :] + dist[i]

            ind = tuple([ndist < dist])
            dist[ind] = ndist[ind]
            prev[ind] = i

        dist[start] = 0
        while any(unvisited):
            unv_ix = np.where(unvisited)[0]
            i = unv_ix[dist[unv_ix].argmin()]
            visit(i)
            unvisited[i] = False

        return prev

    def init_paths(self, first_node, last_node):
        self.first_node = first_node
        self.last_node = last_node
        self.before = self._get_paths(self.first_node, self.w)
        self.after = self._get_paths(self.last_node, self.w.transpose())
        # Load up and remove edges from the ref path

        #path = [self.first_node]
        cur = self.first_node
        last_cur = None
        while self.after[cur] != -1:
            cur = self.after[cur]
            #path.append(cur)
            if last_cur:
                self.edge_set.remove((last_cur, cur))
                log.debug("Removing (%d, %d)", last_cur, cur)
            last_cur = cur

        #self.ref_path = path

    def _get_shortest(self, a, b):
        # Returns the shortest path passing through edge (a, b)

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
        log.debug("%d edges in non-ref edge set.", len(self.edge_set))

        all_paths = set()
        for (i, j) in self.edge_set:
            log.debug("Computing shortest path through edge: (%d, %d)", i, j)
            path = self._get_shortest(i, j)
            if path:
                all_paths.add(tuple(path))

        return list(all_paths)
