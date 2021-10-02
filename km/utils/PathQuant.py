#                             -*- Mode: Python -*-
# PathQuant.py
#

import numpy as np
import logging as log
import sys


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

    @staticmethod
    def output_header():
        sys.stdout.write("\t".join([
                "Database",
                "Query",
                "Type",
                "Variant_name",
                "rVAF",
                "Expression",
                "Min_coverage",
                "Start_offset",
                "Sequence",
                "Reference_expression",
                "Reference_sequence",
                "Info"
            ]) + "\n"
        )


class PathQuant:
    def __init__(self, all_paths, counts):
        self.all_paths = all_paths
        self.nb_kmer = len(counts)
        self.nb_seq = len(all_paths)

        self.counts = np.array(counts, dtype=np.float32)
        self.contrib = np.zeros((self.nb_kmer, self.nb_seq), dtype=np.int32)
        for seq_i, seq in enumerate(all_paths):
            for i in seq:
                self.contrib[i, seq_i] += 1
            # note: self.contrib[seq, seq_i] += 1 would not work in the case of ITDs

        self.coef = None
        self.rVAF = None

        log.debug("%d sequence(s) are observed.", self.nb_seq)

    def compute_coef(self):
        # Set coefficient to zero if all paths use a kmer with 0 coverage
        #if self.contrib[self.counts == 0].min(axis=1).max() == 1:
        #    self.coef = np.zeros((len(c), 1), dtype=np.float32)
        #    return
        (coef, residual, rank, s) = np.linalg.lstsq(self.contrib, self.counts, rcond=None)
        self.coef = coef
        log.debug("Linear fitting = %s", self.coef)

    def refine_coef(self):
        # if max(self.coef) == 0: return
        # applies a gradient descent to get rid of negative coefficients
        self.coef[self.coef < 0] = 0
        last_max_grad = np.inf

        # convergence threshold
        num_iter = 0
        while last_max_grad > 0.01:
            counts_hat = np.dot(self.contrib, self.coef)
            grad = 2 * (self.counts - counts_hat) * self.contrib.T
            grad = grad.sum(axis=1) / self.nb_kmer
            self.coef += 0.1 * grad
            grad[self.coef < 0] = 0
            self.coef[self.coef < 0] = 0
            last_max_grad = np.max(np.abs(grad))
            num_iter += 1
            log.debug(
                "Iteration = %d, max_gradient = %f",
                num_iter,
                last_max_grad
            )
        log.debug("Refined fitting = %s", self.coef)

    def get_ratio(self):
        if max(self.coef) == 0:
            self.rVAF = self.coef
        else:
            self.rVAF = self.coef / np.sum(self.coef)
        return self.rVAF

    def adjust_for_reference(self):
        self.rVAF[0] = np.nan
        self.rVAF[1] = np.nan
        self.coef[self.coef >= 0] = min(self.counts)

    def output(self, db_f, ref_name, name_f, seq_f):
        for i in range(self.nb_seq):
            # if self.rVAF[i] > 0:
            print("%s\t%s\t%s\t%.3f\t%.1f\t%s" % (db_f,
                                                  ref_name,
                                                  name_f(self.all_paths[i]),
                                                  self.rVAF[i], self.coef[i],
                                                  seq_f(self.all_paths[i])))

    def get_paths(self, db_f, ref_name, name_f, seq_f, ref_path, info="",
                  get_min_f=lambda path: 0, start_off=0):
        paths = []
        ref_i = -1
        for i in range(self.nb_seq):
            if list(self.all_paths[i]) == list(ref_path):
                ref_i = i
        for i in range(self.nb_seq):
            if i != ref_i:
                p = Path(db_f, ref_name,
                         name_f(self.all_paths[i]),
                         self.rVAF[i], self.coef[i],
                         get_min_f(self.all_paths[i]),
                         start_off,
                         seq_f(self.all_paths[i]),
                         self.rVAF[ref_i], self.coef[ref_i],
                         seq_f(self.all_paths[ref_i]), info)
                paths += [p]
        return paths
