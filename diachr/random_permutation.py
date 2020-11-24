from collections import defaultdict
from scipy.stats import binom
from collections import defaultdict
from statistics import mean, stdev
import multiprocessing as mp
import itertools
import os

from .diachr_util import calculate_binomial_p_value, find_minimum_powered_n
from .diachromatic_parser import DiachromaticParser
from .diachromatic_interaction import DiachromaticInteraction

class RandomPermutation:
    def __init__(self, combined_interaction_file: str, alpha: float, threads: int = 1, prefix:str="PERMUTATION") -> None:
        self._pval_memo = defaultdict(float)
        self._nominal_alpha = alpha
        self._N_DICT = defaultdict(int)
        self._nthreads = threads
        self._prefix = prefix
        if not os.path.isfile(combined_interaction_file):
            raise ValueError("Need to provide path to a single combined interaction file")
        parser = DiachromaticParser(combined_interaction_file)
        self._interaction_d = parser.get_interaction_dict()
        for _, interaction in self._interaction_d.items():
            self._N_DICT[interaction.total_readpairs] += 1
        print("[INFO] We ingested {} interactions from {}".format(len(self._interaction_d), combined_interaction_file))
        # Minimum number of read pairs required for significance
        self._minimum_powered_n, self._pv_of_minimum_powered_n = find_minimum_powered_n(self._nominal_alpha)
        # Total number of input interactions
        self._n_interaction = 0
        # Number of interaction that cannot be significant
        self._n_underpowered_interaction = 0
        self._n_directed_interaction = 0
        self._n_undirected_interaction = 0
        # Calculate the distribution of p-values for the observed data
        for _, interaction in self._interaction_d.items():
            self._n_interaction += 1
            # Skip interactions that cannot be significant
            n = interaction.total_readpairs
            if n < self._minimum_powered_n:
                self._n_underpowered_interaction += 1
                continue
            key = interaction.key
            if key in self._pval_memo:
                pv = self._pval_memo[key]
            else:
                pv = calculate_binomial_p_value(interaction.simple, interaction.twisted)
                self._pval_memo[key] = pv
            # Count interaction as directed or undirected
            if pv < self._nominal_alpha:
                self._n_directed_interaction += 1
            else:
                self._n_undirected_interaction += 1
            self._N_DICT[interaction.total_readpairs] += 1
        # Output some statistics about original interactions
        print("[INFO] Total number interactions: {}".format(len(self._interaction_d)))
        print("\t[INFO] Number of underpowered interactions (discarded): {}".format(self._n_underpowered_interaction))
        print("\t[INFO] Number of undirected interactions: {}".format(self._n_undirected_interaction))
        print("\t[INFO] Number of interactions with nominally significant P-values (<{}): {}".format(self._nominal_alpha, self._n_directed_interaction))



    def _binomial_p_value(self, simple_count: int, twisted_count: int) -> float:
        """
        Calculation of the binomial P-value for deviation of simple/twisted counts from expected p=0.5
        :param simple_count: Number of simple read pairs
        :param twisted_count: Number of twisted read pairs
        :return: Binomial P-value
        """
        # Create key from simple and twisted read pair counts
        key = "{}-{}".format(simple_count, twisted_count)
        if key in self._pval_memo:
            return self._pval_memo[key]
        else:
            p_value = calculate_binomial_p_value(simple_count, twisted_count)
            self._pval_memo[key] = p_value
            return p_value

    def _perform_one_iteration(self):
        """
        This function performs one iteration of the random permutation analysis. It permutes the simple and twisted read
        pair counts of all interactions and then determines and returns the number of significant interactions.

        :return: Number of significant interactions after permutation
        """

        # Init number of significant interactions after permutation
        w_perm_significant = 0

        # Iterate dictionary with numbers of interactions foreach read pair number n
        for n, i_num in self._N_DICT.items():
            # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p = 0.5, size = i_num))

            for simple_count in simple_count_list:
                twisted_count = n - simple_count
                key = "{}-{}".format(simple_count, twisted_count)
                if key in self._pval_memo:
                    pv = self._pval_memo[key]
                else:
                    pv = calculate_binomial_p_value(simple_count, twisted_count)
                    self._pval_memo[key] = pv
                # Count significant interactions
                if pv < self._nominal_alpha:
                    w_perm_significant += 1

        return w_perm_significant

    def _perform_m_iterations(self, m_iter: int):
        w_sig_p_list = []
        print("\t[INFO] Batch: Performing " + str(m_iter) + " iterations ...")
        for i in range(1, m_iter + 1):
            w_sig_p_list.append(self._perform_one_iteration())
        return w_sig_p_list

    def permute(self, ITER_NUM: int) -> None:
        # Perform permutation analysis
        print("[INFO] Performing permutation analysis with " + str(ITER_NUM) + " iterations ...")
        pool = mp.Pool(processes = self._nthreads )

        # Perform permutation for 'thread_num' batches with balanced numbers of iterations
        batch_iter_nums = list(itertools.repeat(int(ITER_NUM / self._nthreads), self._nthreads))
        if 0 < ITER_NUM - self._nthreads * int(ITER_NUM / self._nthreads):
            # If there is a remainder, add it to the first element of list
            batch_iter_nums[0] = batch_iter_nums[0] + ITER_NUM - self._nthreads * int(ITER_NUM / self._nthreads)
        results = [pool.apply_async(self._perform_m_iterations, args=(batch_iter_num,)) for batch_iter_num in batch_iter_nums]
        pool.close()

        # Combine results from different batches
        batch_results = [p.get() for p in results]
        print("[INFO] Combining results from different batches ...")

        # List with numbers of significant interactions after permutation for each iteration
        w_sig_p_list = list(itertools.chain.from_iterable(batch_results))
        random_better_than_observed = len([i for i in w_sig_p_list if i > self._n_directed_interaction])
        # Calculate average number of significant permuted interactions
        w_sig_p_average = mean(w_sig_p_list)
        # Calculate standard deviation of significant permuted interactions
        w_sig_p_sd = stdev(w_sig_p_list)
        z_score = "{0:.2f}".format((self._n_directed_interaction - w_sig_p_average) / w_sig_p_sd)
        total = self._n_directed_interaction + self._n_undirected_interaction
        percentage_observed = "{0:.2f}".format(self._n_directed_interaction / total)
        percentage_permuted = "{0:.2f}".format(w_sig_p_average / total)
        print("\t[INFO] OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_INDEFINABLE_INTERACTION"
               "\tN_UNDIRECTED_INTERACTION"
               "\tN_DIRECTED_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_INTERACTION"
               "\tSD_PERMUTED_DIRECTED_INTERACTION"
               "\tZ_SCORE"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED")

        print("\t[INFO] " + self._prefix + "\t"
            + str(ITER_NUM) + "\t"
            + str(self._nominal_alpha) + "\t"
            + str(self._minimum_powered_n) + "\t"
            + str(self._n_interaction) + "\t"
            + str(self._n_underpowered_interaction) + "\t"
            + str(self._n_undirected_interaction) + "\t"
            + str(self._n_directed_interaction) + "\t"
            + "{0:.2f}".format(w_sig_p_average) + "\t"
            + "{0:.2f}".format(w_sig_p_sd) + "\t"
            + str(z_score) + "\t"
            + str(random_better_than_observed) + "\n")

        print("{} out of {} permutations had more signficant p values than in the observed data.".format(random_better_than_observed, str(ITER_NUM)))
        # Write results to the file
        # Write table with parameters and summary statistics to file
        # Prepare stream for output of summary statistics
        tab_file_stats_output = self._prefix + "_permutation_stats.txt"
        fh = open(tab_file_stats_output, 'wt')
        fh.write("OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_INDEFINABLE_INTERACTION"
               "\tN_UNDIRECTED_INTERACTION"
               "\tN_DIRECTED_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_INTERACTION"
               "\tSD_PERMUTED_DIRECTED_INTERACTION"
               "\tZ_SCORE"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED\n")

        fh.write("\t[INFO] " + self._prefix + "\t"
            + str(ITER_NUM) + "\t"
            + str(self._nominal_alpha) + "\t"
            + str(self._minimum_powered_n) + "\t"
            + str(self._n_interaction) + "\t"
            + str(self._n_underpowered_interaction) + "\t"
            + str(self._n_undirected_interaction) + "\t"
            + str(self._n_directed_interaction) + "\t"
            + "{0:.2f}".format(w_sig_p_average) + "\t"
            + "{0:.2f}".format(w_sig_p_sd) + "\t"
            + str(z_score) + "\t"
            + str(random_better_than_observed) + "\n")

        fh.close()

        # Write numbers of significant interactions for all iterations to file
        # Prepare stream for output of significant interaction numbers for all iterations (w_i)
        txt_file_w_sig_p_output = self._prefix + "_permutation_w_sig_p.txt"
        fh = open(txt_file_w_sig_p_output, 'wt')
        for w_sig_p in w_sig_p_list:
            fh.write(str(w_sig_p) + "\n")
        fh.close()