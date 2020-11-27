#                             -*- Mode: Python -*-
# PathQuant.py
#

import numpy as np
import logging as log


class Path:
    def __init__(self, db_f, ref_name, variant_name, ratio, expression,
                 min_coverage, start_off, sequence, ref_ratio, ref_expression,
                 ref_sequence, note):
        self.db_name = db_f
        self.ref_name = ref_name
        self.variant_name = variant_name
        self.rVAF = ratio
        self.expression = expression
        self.min_coverage = min_coverage
        self.start_off = start_off
        self.sequence = sequence
        self.ref_ratio = ref_ratio
        self.ref_expression = ref_expression
        self.ref_sequence = ref_sequence
        self.note = note

    def __str__(self):
        return "%s\t%s\t%s\t%.3f\t%.1f\t%d\t%d\t%s\t%.1f\t%s\t%s" % (
            self.db_name,
            self.ref_name,
            self.variant_name,
            self.rVAF,
            self.expression,
            self.min_coverage,
            self.start_off,
            self.sequence,
            self.ref_expression,
            self.ref_sequence,
            self.note)

    def __list__(self):
        return self.__str__().split('\t')

    def __getitem__(self, i):
        return self.__list__()[i]

    def get_min_cov(self):
        return self.min_coverage

    def get_sequence(self):
        return self.sequence

    def get_variant_name(self):
        return self.variant_name


class PathQuant:
    def __init__(self, all_path, counts):
        self.all_path = all_path
        self.nb_kmer = len(counts)
        self.nb_seq = len(all_path)

        self.contrib = np.zeros((self.nb_kmer, self.nb_seq), dtype=np.int32)
        self.counts = np.zeros((self.nb_kmer, 1), dtype=np.float32)

        self.coef = None
        self.rVAF = None

        seq_i = 0
        log.debug("%d sequence(s) are observed.", self.nb_seq)

        for s in all_path:
            for i in s:
                self.contrib[i, seq_i] += 1
            seq_i += 1

        self.counts[:, 0] = counts

    def compute_coef(self):
        # Set coefficient to zero if all paths use a kmer with 0 coverage
        # for c,d in zip(self.contrib,self.counts):
        #     print c, d
        #     if min(c) == 1 and d == 0:
        #         self.coef = np.zeros((len(c), 1), dtype=np.float32)
        #         return
        (coef, residual, rank, s) = np.linalg.lstsq(self.contrib, self.counts, rcond=None)
        self.coef = coef
        log.debug("Linear fitting = %s", self.coef.flatten())

    def refine_coef(self):
        # if max(self.coef) == 0: return
        # applies a gradient descent to get rid of negative coefficients
        self.coef[self.coef < 0] = 0
        last_max_grad = np.inf
        num_iter = 0

        # convergence threshold
        while last_max_grad > 0.01:
            grad = np.zeros_like(self.coef, dtype=np.float32)
            counts_hat = np.dot(self.contrib, self.coef)
            for j in range(self.nb_seq):
                grad[j, 0] = np.sum(2 * (self.counts - counts_hat) *
                                    self.contrib[:, j].reshape(counts_hat.shape))
            grad /= self.nb_kmer
            self.coef += 0.1 * grad
            grad[self.coef < 0] = 0
            self.coef[self.coef < 0] = 0
            last_max_grad = np.max(np.abs(grad))
            num_iter += 1
            log.debug("Iteration = %d, max_gradient = %f",
                      num_iter,
                      last_max_grad)
        log.debug("Refined fitting = %s", self.coef.flatten())

    def get_ratio(self):
        if max(self.coef) == 0:
            self.rVAF = self.coef
        else:
            self.rVAF = self.coef / np.sum(self.coef)
        return self.rVAF

    def adjust_for_reference(self):
        self.rVAF[0] = np.nan
        self.rVAF[1] = np.nan
        self.coef[self.coef >= 0] = min(self.counts)[0]

    @staticmethod
    def output_header():
        print("Database\tQuery\tType\tVariant_name\trVAF\tExpression\tMin_coverage\tStart_offset\tSequence\tReference_expression\tReference_sequence\tInfo")

    def output(self, db_f, ref_name, name_f, seq_f):
        for i in range(self.nb_seq):
            # if self.rVAF[i] > 0:
            print("%s\t%s\t%s\t%.3f\t%.1f\t%s" % (db_f,
                                                  ref_name,
                                                  name_f(self.all_path[i]),
                                                  self.rVAF[i], self.coef[i],
                                                  seq_f(self.all_path[i])))

    def get_paths(self, db_f, ref_name, name_f, seq_f, ref_path, info="",
                  get_min_f=lambda path: 0, start_off=0):
        paths = []
        ref_i = -1
        for i in range(self.nb_seq):
            if list(self.all_path[i]) == list(ref_path):
                ref_i = i
        for i in range(self.nb_seq):
            if i != ref_i:
                p = Path(db_f, ref_name,
                         name_f(self.all_path[i]),
                         self.rVAF[i], self.coef[i],
                         get_min_f(self.all_path[i]),
                         start_off,
                         seq_f(self.all_path[i]),
                         self.rVAF[ref_i], self.coef[ref_i],
                         seq_f(self.all_path[ref_i]), info)
                paths += [p]
        return paths
