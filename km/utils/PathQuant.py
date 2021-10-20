#                             -*- Mode: Python -*-
# PathQuant.py
#

import numpy as np
import logging as log
import sys


class Path:
    """Path object. Wraps path data into a printable format.

    Methods
    -------
    get_min_cov
    get_sequence
    get_variant_name
    output_header
    """

    def __init__(self, db_f, ref_name, variant_name, rvaf,
                 expression, ref_expression, min_coverage,
                 ref_min_coverage, sequence, ref_sequence, note):
        self.db_name = db_f
        self.ref_name = ref_name
        self.variant_name = variant_name
        self.rvaf = rvaf
        self.expression = expression
        self.ref_expression = ref_expression
        self.min_coverage = min_coverage
        self.ref_min_coverage = ref_min_coverage
        self.sequence = sequence
        self.ref_sequence = ref_sequence
        self.note = note

    def __str__(self):
        return "%s\t%s\t%s\t%.3f\t%.1f\t%.1f\t%d\t%d\t%s\t%s\t%s" % (
            self.db_name,
            self.ref_name,
            self.variant_name,
            self.rvaf,
            self.expression,
            self.ref_expression,
            self.min_coverage,
            self.ref_min_coverage,
            self.sequence,
            self.ref_sequence,
            self.note)

    def __list__(self):
        return self.__str__().split('\t')

    def __getitem__(self, i):
        return self.__list__()[i]

    def get_min_cov(self):
        """Return minimum coverage"""

        return self.min_coverage

    def get_sequence(self):
        """Return variant sequence"""

        return self.sequence

    def get_variant_name(self):
        """Return variant description (also called variant name)"""

        return self.variant_name

    @staticmethod
    def output_header():
        """Output pre-formatted header"""

        sys.stdout.write("\t".join([
                "Database",
                "Query",
                "Type",
                "Variant",
                "rVAF",
                "Expression",
                "Ref_expression",
                "Min_coverage",
                "Ref_min_coverage",
                "Sequence",
                "Reference_sequence",
                "Info"
            ]) + "\n"
        )


class PathQuant:
    """Paths from graph analysis are quantified according
    to their respecting kmer counts.

    An expectation-maximization (EM) function is used to
    iteratively compute (E step) and maximize (M step) the
    likelihood function. The EM algorithm will always converge
    towards the (global) maximum likelihood (ML) value. This is
    due to the fact that ML optimization has only one local maximum.
    Note that we use the log-likelihood as a proxy for the
    likelihood computation as maximizing either one is equivalent
    to maximizing the other, with the added benefit of having more
    precise float values with the log-likelihood.

    The EM will compute a ratio for the input paths.

    Note
    ----
    Disclaimer: In the context of kmer-based quantification, the
    ML maximization problem is limited by the fact that kmers are
    not independently sampled. For example, with k=31, a read
    of 100bp will generate 70 kmers. Furthermore, sequence lengths
    for ML estimation are limited to the kmer walking paths
    and completely dependent on the original target sequence. Due
    to these concerns, the quantification shown in km should be used
    carefully and treated as an estimation at best. Significant
    conclusions should not be drawn from these values.

    Attributes
    ----------
    all_paths : (N,) list
        A list of all paths that need quantification.
    counts : (M,) array_like
        A 1-dimensional vector with kmer counts.
    contrib : (M, N) array_like
        A 2-dimension array with 1s and 0s to indicate
        which paths contain which kmer.
    nb_seq : int
        Number of path sequences to quantify (N).
    nb_kmer : int
        Number of kmers used for quantification (M).
    ratio : (N,) array_like
        Ratio values for all paths. Always sums to 1.

    Methods
    -------
    quantify
    """

    def __init__(self, all_paths, counts):
        self.all_paths = all_paths
        self.counts = np.array(counts, dtype=np.float32)

        self.nb_seq = len(all_paths)
        self.nb_kmer = len(counts)

        self.contrib = np.zeros((self.nb_kmer, self.nb_seq), dtype=np.int32)
        for seq_i, seq in enumerate(all_paths):
            for i in seq:
                self.contrib[i, seq_i] += 1
            # note: self.contrib[seq, seq_i] += 1 would not work in the case of ITDs

        self.ratio = None
        self.expression = None

        log.info("%d sequence(s) are observed.", self.nb_seq)

    def quantify(self):
        """Perform EM quantification.

        Variable `rho` represents relative abundances of transcripts.
        Variable `alpha` represents relative kmer allocation to each
        transcript/path. When lengths are equal, rho == alpha. Of
        note, `rho` can be deduced from alpha and sequence lengths
        (and vice-versa). Grossly speaking, I find it useful to
        think of `rho` as the TPM value and `alpha` as read counts.

        Note that computing the log-likelihood does nothing other
        than asserting that the EM actually converges. Chances are
        that if we remove its computation, this function will
        continue to behave exactly the same, but it is not very
        expensive to compute and gives an additional assurance for
        the method.

        Note
        ----
        By virtue of its nature, EM converges towards the global ML
        maximum, so we could start off with any rho (ratio) values.
        However, in the case of 2 identical sequences being compared
        (when quantifying the reference, for instance), ratio values
        will remain the same as the starting values. By using equal
        ratio values for all sequences we completely avoid this
        problem.

        See also
        -------
        https://data-science-sequencing.github.io/Win2018/lectures/lecture12/
        """

        nn_ix = np.where(self.contrib.sum(axis=1) > 0)[0]
        all_counts = self.counts[nn_ix].sum()
        paths_lengths = self.contrib.sum(axis=0)
        # no need to compute - 31 + 1 because we are counting kmers not nucleotides

        def calc_em(cur_rho):
            # compute expression
            w = self.contrib[nn_ix,:] * cur_rho
            weighted_contrib = w.T / np.maximum(w.sum(axis=1), np.finfo(float).eps)
            #^ divide non-zero values by sum and all zero values by eps to return 0
            counts = weighted_contrib * self.counts[nn_ix]
            count_sums = counts.sum(axis=1)
            expression = count_sums / paths_lengths

            with np.errstate(divide='raise', invalid='raise'):
                assert np.isclose(all_counts, count_sums.sum())
                try:
                    # compute alpha
                    alpha = count_sums / all_counts
                    # compute rho
                    new_rho = alpha / paths_lengths
                    new_rho = new_rho / new_rho.sum()  # note that new_rho /= new_rho.sum() would
                except FloatingPointError:             # get applied before exception catching
                    assert count_sums.sum() == 0
                    alpha = new_rho = count_sums
                    #log.info('Counts are all 0.')

            # compute log-likelihood
            selection_prob = np.dot(self.contrib[nn_ix], alpha / paths_lengths)
            counts = self.counts[nn_ix]
            if len(np.where(selection_prob == 0)[0]):
                nonzero = np.where(selection_prob > 0)[0]
                selection_prob = selection_prob[nonzero]
                counts = counts[nonzero]
            log_likelihood = np.sum(np.log(selection_prob) * counts)
            return new_rho, alpha, expression, log_likelihood

        rho = np.array([1/self.nb_seq]*self.nb_seq)

        log_likelihood = -np.inf
        prev_rho = np.array([-np.inf]*rho.shape[0])
        num_iter = 0
        done = False

        while max(abs(rho - prev_rho)) >= 1e-6:
            num_iter += 1
            prev_rho = rho
            prev_log_likelihood = log_likelihood

            rho, alpha, expression, log_likelihood = calc_em(rho)

            try:
                assert log_likelihood >= prev_log_likelihood
                #^ 'or equal' for cases where at least one alpha value is 0
            except AssertionError:
                log.info("ERROR @ Log likelihood assertion failed")

            log.debug(
                "Iteration=%d, rho=%s, LogLikelihood=%f",
                num_iter,
                str(rho.round(3)),
                log_likelihood
            )

        log.info("EM quantification = %s", rho)

        self.ratio = rho
        self.expression = expression

        #w = self.contrib[nn_ix,:] * rho
        #weighted_contrib = w.T / np.maximum(w.sum(axis=1), np.finfo(float).eps)
        #counts = weighted_contrib * self.counts[nn_ix]
        #log_counts = np.log(np.where(counts == 0, 1, counts))
        #self.geomexpr = np.exp(log_counts.sum(axis=1) / paths_lengths)
        ## special case where all counts are 0: geometric mean needs adjustment
        #self.geomexpr[counts.sum(axis=1) == 0] = 0


    def adjust_for_reference(self):
        self.ratio = np.array([1, 0])
        self.expression = np.array([self.expression.sum(), 0])


    #def compute_coef(self):
    #    # Set coefficient to zero if all paths use a kmer with 0 coverage
    #    #if self.contrib[self.counts == 0].min(axis=1).max() == 1:
    #    #    self.coef = np.zeros((len(c), 1), dtype=np.float32)
    #    #    return
    #    (coef, residual, rank, s) = np.linalg.lstsq(self.contrib, self.counts, rcond=None)
    #    self.coef = coef
    #    log.debug("Linear fitting = %s", self.coef)

    #def refine_coef(self):
    #    # if max(self.coef) == 0: return
    #    # applies a gradient descent to get rid of negative coefficients
    #    self.coef[self.coef < 0] = 0
    #    last_max_grad = np.inf

    #    # convergence threshold
    #    num_iter = 0
    #    while last_max_grad > 0.01:
    #        counts_hat = np.dot(self.contrib, self.coef)
    #        grad = 2 * (self.counts - counts_hat) * self.contrib.T
    #        grad = grad.sum(axis=1) / self.nb_kmer
    #        self.coef += 0.1 * grad
    #        grad[self.coef < 0] = 0
    #        self.coef[self.coef < 0] = 0
    #        last_max_grad = np.max(np.abs(grad))
    #        num_iter += 1
    #        log.debug(
    #            "Iteration = %d, max_gradient = %f",
    #            num_iter,
    #            last_max_grad
    #        )
    #    log.debug("Refined fitting = %s", self.coef)

    #def get_ratio(self):
    #    if max(self.coef) == 0:
    #        self.ratio = self.coef
    #    else:
    #        self.ratio = self.coef / np.sum(self.coef)
    #    return self.ratio

    #def adjust_for_reference(self):
    #    self.ratio[0] = np.nan
    #    self.ratio[1] = np.nan
    #    self.coef[self.coef >= 0] = min(self.counts)
