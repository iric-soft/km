#                             -*- Mode: Python -*-
# Graph.py --- Generic graph class (directed, weighted)
#

import logging as log
import numpy as np
import time
import datetime

import resource, sys
resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
sys.setrecursionlimit(10**6)

class Graph:
    # Nodes are identified by 0..(n-1)

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

    def __getitem__(self, (i, j)):
        return self.w[i, j]

    def __setitem__(self, (i, j), new_val):
        self.w[i, j] = new_val
        self.edge_set.add((i, j))

    def _get_paths(self, start, w):
        prev = np.empty(self.n, dtype=np.int32)  # represents indexes of previous kmer
        prev.fill(-1)
        dist = np.empty(self.n, dtype=np.float32)  # represents passed over kmers with their order
        dist.fill(np.inf)
        unvisited = set(range(self.n))  # represents index of kmers that weren't passed over yet
        
        def visit(i):
            # i is previous; j is current
            # all nodes in graph that continue i, with their weights + weight of i
            ndist = w[i, :] + dist[i]
            # iterate through all nodes
            memory = 0
            for j in range(self.n):
                # find the kmers that continue i (ignoring already passed over values (from + dist[i]))
                if ndist[j] < dist[j]:  # will only work when dist[j] == inf
                    if memory:
                        log.debug("Branching: kmers=(%d,%d); weights=(%.3f,%.3f,%.3f)"
                                  % (i, j, memory, ndist[j], dist[j]))
                    dist[j] = ndist[j]
                    memory = ndist[j]
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
        
        # This loop is the bottleneck in find_fusion (in addition to path quantification in some cases)
        # Finds shortest path by keeping path counts and order in the dist list
        #time_start = time.time()
        while unvisited:
            i = min_unvisited()  # should be the last kmer that was visited
            unvisited.remove(i)
            visit(i)
        #print str(datetime.timedelta(seconds=time.time() - time_start))
        
        return prev
        

    def init_paths(self, first_node, last_node):
        # Updates self.before, self.after and self.edge_set before running all_shortest
        self.first_node = first_node
        self.last_node = last_node
        log.debug("Initialize path (before)")
        self.before_shortest = self._get_paths(self.first_node, self.w)  # contains idx of previous kmer
        log.debug("Initialize path (after)")
        self.after_shortest = self._get_paths(self.last_node, self.w.transpose())  # idx of next kmer
        
        self.before = np.empty(self.n, dtype=np.int32)
        self.after = np.empty(self.n, dtype=np.int32)
        self.after_ref = np.empty(self.n, dtype=tuple)
        
        for i in range(self.n):
            bef = set(np.where(self.after_shortest == i)[0])
            bef.add(self.before_shortest[i])
            if self.first_node in bef:
                self.before[i] = self.first_node
            else:
                for b in list(bef):
                    ref_branch = 0
                    if np.isclose(self.w[b, i], 0.001):
                        ref_branch += 1
                if ref_branch <= 1:  # Reference is unambiguous, trust shortest paths
                    self.before[i] = self.before_shortest[i]
                else:
                    self.before[i] = -1
            
            aft = set(np.where(self.before_shortest == i)[0])
            aft.add(self.after_shortest[i])
            if self.last_node in aft:
                self.after[i] = self.last_node
            else:
                for a in list(aft):
                    ref_branch = 0
                    if np.isclose(self.w[i, a], 0.001):
                        ref_branch += 1
                if ref_branch <= 1:  # Reference is unambiguous, trust shortest paths
                    self.after[i] = self.after_shortest[i]
                else:
                    self.after[i] = -1
            
            aft = set(np.where(self.before_shortest == i)[0])
            aft.add(self.after_shortest[i])  # When the next kmer has many previous ones
            for a in list(aft):
                if self.w[i, a] > 0.99:
                    aft.remove(a)
            self.after_ref[i] = tuple(aft)

        # Load up and remove edges from the ref path
        # Consider implementing a tree instead of a recursive function
        def remove_edges(first_local_node):
            cur = first_local_node
            last_cur = cur
            while self.after_ref[cur] and self.after_ref[cur][0] != -1 and cur not in visited:
                visited.add(cur)
                cur = self.after_ref[cur]
                for c in cur:
                    if (last_cur, c) in self.edge_set:
                        self.edge_set.remove((last_cur, c))
                        log.debug("Removing(%d, %d)", last_cur, c)
                if len(cur) > 1:
                    log.debug("Branching ...")
                    for c in cur:
                        remove_edges(c)
                    break
                else:
                    cur = cur[0]
                    last_cur = cur

        visited = set()
        remove_edges(self.first_node)
    
    
    def get_path_score(self, path):
        # Returns the shortest path passing through edge (a, b)
        score = 0
        for i in range(len(path)-1):
            score += self.w[path[i], path[i+1]]
        
        return score

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
        # Only keep paths from source to sink
        if path[0] != self.first_node or path[-1] != self.last_node:
            return None
        
        return(path)

    def all_shortest(self) :
        all_paths = set()
        log.debug("%d edges in non-ref edge set.", len(self.edge_set))
        for (i, j) in self.edge_set:
            log.debug("Computing shortest path through edge: (%d, %d)", i, j)
            path = self.get_shortest(i, j)
             
            if path:
                all_paths.add(tuple(path))
       
        return list(all_paths)

    def diff_path_without_overlap(self, ref, seq, k):
        # Returns:
        #  - start
        #  - stop_ref
        #  - stop_variant
        #  - kmers_ref
        #  - kmers_variant
        #  - stop_ref_fully_trimmed
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

        return(i, j_ref, j_seq, ref[i:j_ref], seq[i:j_seq], k_ref)
