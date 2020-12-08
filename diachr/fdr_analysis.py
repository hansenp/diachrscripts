from collections import defaultdict
from scipy.stats import binom
import numpy as np
import math
import gzip

from .diachr_util import calculate_binomial_logsf_p_value
from diachr.diachromatic_interaction import DiachromaticInteraction
from diachr.diachromatic_interaction_parser import DiachromaticInteractionParser

class FdrAnalysis:
    """
    Initially, a Diachromatic interaction file is ingested and a binomial P-value is calculated for each interaction and
    stored in a list. Furthermore, the total numbers of interactions for given read pair numbers n (n is the sum of simple
    and twisted read pairs) are stored in dictionary with n as keys and interaction numbers as values.

    This dictionary is then used to generate a list of random P-values for all interactions. A random P-value for an
    interaction with n read pairs is generated by drawing a number of simple read pairs s' from a binomial distribution
    B(n, p = 0.5), setting the number of twisted read pairs to t' to n - s' and calculating the corresponding binomial
    P-value.

    Finally, the FDR is estimated for increasing P-value thresholds. For each threshold, the number of significant
    interactions is determined for the list of original P-values (S_o) and for the lists of randomized P-values (S_p) and
    S_p / S_o is used as estimator for the FDR.
    """

    def __init__(self, diachromatic_interaction_file:str, 
                        fdr_threshold:float=0.25, 
                        p_val_c_min:float=0.00025, 
                        p_val_c_max:float=0.05, 
                        step_size:float=0.00025, 
                        prefix:str='OUTPREFIX') -> None:
        self._inputfile = diachromatic_interaction_file
        if fdr_threshold <= 0 or fdr_threshold > 1.0:
            raise ValueError("FDR threshold must be in (0,1)")
        self._fdr_threshold = fdr_threshold
        if p_val_c_min <= 0 or p_val_c_min > 1.0:
            raise ValueError("p_val_c_min must be in (0,1)")
        self._p_val_c_min = p_val_c_min
        if p_val_c_max <= 0 or p_val_c_max > 1.0:
            raise ValueError("p_val_c_max must be in (0,1)")
        if p_val_c_max < p_val_c_min:
            raise ValueError("p_val_c_max must be larger than p_val_c_min")
        self._p_val_c_max = p_val_c_max
        self._step_size = step_size
        self._prefix = prefix
        self._pval_memo = defaultdict(float)
        print("[INFO] " + "Input parameters")
        print("\t[INFO] --out_prefix: " + self._prefix)
        print("\t[INFO] --diachromatic-interaction-file: " + self._inputfile)
        print("\t[INFO] --fdr-threshold: " + str(self._fdr_threshold))
        print("\t[INFO] --p-val-c-min: " + str(self._p_val_c_min))
        print("\t[INFO] --p-val-c-max: " + str(self._p_val_c_max))
        print("\t[INFO] --p-val-step-size: " + str(self._step_size))
        # Dictionary that stores the numbers of interactions with n read pairs
        self._N_DICT = defaultdict(int)
        # Total number of input interactions
        self._n_interaction = 0
        # List of P-values for observed interactions
        self._p_val_o_list = []

        # Get list of Dichromatic interaction objects
        parser = DiachromaticInteractionParser()
        parser.parse_file(diachromatic_interaction_file)
        d_inter_list = parser.i_list
        np.random.seed(42)
        n_progress = 0
        for d_inter in d_inter_list:
            n_progress += 1
            if n_progress % 100000 == 0:
                print("\t[INFO] Processed " + str(n_progress) + " interaction objects ...")

            # Get neg_log_p_value
            neg_log_p_value = float("{:.2f}".format(-calculate_binomial_logsf_p_value(d_inter.n_simple, d_inter.n_twisted)))
            self._p_val_o_list.append(neg_log_p_value)

            # Add the sum of simple and twisted read pair counts to dictionary that will be used for randomization
            if d_inter.rp_total in self._N_DICT:
                self._N_DICT[d_inter.rp_total] += 1
            else:
                self._N_DICT[d_inter.rp_total] = 1

        self._n_interaction = n_progress
        print("[INFO] Total number of interactions: {}".format(self._n_interaction))
        print("[INFO] Getting sorted lists of P-values ...")

    def _get_pvals_permuted_counts(self):
        """
        This function generates and returns a list of randomized binomial P-values.

        :return: List of randomized P-values
        """
        pvals_permuted_counts = []
        # Iterate dictionary with numbers of interactions foreach read pair number n
        for n, i_num in self._N_DICT.items():
            # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p = 0.5, size = i_num))
            for simple_count in simple_count_list:
                twisted_count = n - simple_count
                # Get binomial P-value
                key = "{}-{}".format(simple_count, twisted_count)
                if key in self._pval_memo:
                    pv_nnl = self._pval_memo[key]
                else:
                    pv_nnl = self._binomial_nlogsf_p_value(simple_count, twisted_count)
                    self._pval_memo[key] = pv_nnl
                pvals_permuted_counts.append(pv_nnl)
        return pvals_permuted_counts


    def _binomial_nlogsf_p_value(self, simple_count:int, twisted_count:int):
        """
        Locally defined method for the calculation of the binomial P-value that uses a dictionary that keeps track of
        P-values that already have been calculated.
        :param simple_count: Number of simple read pairs
        :param twisted_count: Number of twisted read pairs
        :return: Binomial P-value
        """
        # Create key from simple and twisted read pair counts
        key = "{}-{}".format(simple_count, twisted_count)
        # Check whether a P-value for this combination of simple and twisted counts has been calculated already
        if key in self._pval_memo:
            return self._pval_memo[key]
        else:
            # Calculate P-value and add to dictionary
            nnl_p_value = -1*calculate_binomial_logsf_p_value(simple_count, twisted_count)
            self._pval_memo[key] = nnl_p_value
            return nnl_p_value

    def calculate_fdr_threshild(self):
        fdr_last = 1
        S_o_last = 0
        S_r_last = 0
        pc_last = 0
        filename = self._prefix + "_fdr_analysis_results.txt"
        fh = open(filename, 'wt')
        fh.write("OUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O" + "\n")
        # Get list of randomized P-values
        p_val_r_list = self._get_pvals_permuted_counts()
        self._p_val_o_list.sort(key=float, reverse=True)
        p_val_r_list.sort(key=float, reverse=True)
        print("[INFO] Estimating FDR for increasing P-value thresholds ...")
        # Estimate FDR for increasing P-value thresholds from original and randomized P-value lists
        S_o = 0 # start index in sorted list
        S_r = 0 # start index in sorted list
        for pc in np.arange(self._p_val_c_min, self._p_val_c_max, self._step_size):
        # Transform P-value to negative natural logarithm (nnl)
            pc_nnl = - math.log(pc)
            # Get number of significant interactions for original P-values
            #S_o = sum(pc_nnl < pv_nnl for pv_nnl in p_val_o_list) # slow
            while S_o < len(self._p_val_o_list) and self._p_val_o_list[S_o] > pc_nnl:
                S_o += 1

            # Get number of significant interactions for randomized P-values
            #S_r = sum(pc_nnl < pv_nnl for pv_nnl in p_val_r_list) # slow
            while S_r < len(p_val_r_list) and p_val_r_list[S_r] > pc_nnl:
                S_r += 1

            # Estimate FDR
            fdr = S_r / S_o
            fh.write(self._prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(S_r) + "\t" + str(S_o) + "\n")
            print("\t" + self._prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(S_r) + "\t" + str(S_o))

            # Keep track of the largest P-value that satisfies the FDR threshold
            if fdr < self._fdr_threshold:
                fdr_last = fdr
                S_o_last = S_o
                S_r_last = S_r
                pc_last = pc

        fh.close()
        # Print results for the largest P-value that satisfies the FDR threshold to the screen
        print()
        print("\tOUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O")
        print("\t" + self._prefix + "\t" + str(fdr_last) + "\t" + str(pc_last) + "\t" + str(S_r_last) + "\t" + str(S_o_last))
