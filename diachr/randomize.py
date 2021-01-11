from scipy.stats import binom
from numpy import arange
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

class Randomize:

    def __init__(self, random_seed: int = 42, interaction_set: DiachromaticInteractionSet = None):
        self._random_seed = random_seed
        self._interaction_set = interaction_set

    def get_pval_tresh_at_chosen_fdr_tresh(self,
                                           chosen_fdr_thresh: float = 0.05,
                                           pval_thresh_max: float = 0.05,
                                           pval_thresh_step_size: float = 0.00025
                                           ):

        # Check input parameters
        if chosen_fdr_thresh <= 0 or 1.0 < chosen_fdr_thresh:
            raise ValueError("FDR threshold must be in ]0,1]!")
        if pval_thresh_max <= 0 or 1.0 < pval_thresh_max:
            raise ValueError("Maximum P-value threshold must be in ]0,1]!")
        if pval_thresh_max <= pval_thresh_step_size:
            raise ValueError("Step size must be smaller than maximum P-value threshold!")

        # Get list of observed P-values and dictionary that stores the numbers of interactions with n read pairs
        pvals_observed = []
        for d_inter in self._interaction_set.interaction_list:
            pvals_observed.append(d_inter._nln_pval)

        # Get dictionary that stores the numbers of interactions with n read pairs and list of random P-values
        rp_inter_dict = self._get_rp_inter_dict()
        pvals_randomized = self._get_list_of_p_values_from_randomized_data(rp_inter_dict=rp_inter_dict)

        # Determine numbers of significant interactions in observed and randomized data for increasing P-value
        # thresholds and estimate FDR.
        pvals_observed.sort(key=float, reverse=True)
        pvals_randomized.sort(key=float, reverse=True)

        p_thresh_list = arange(pval_thresh_step_size, pval_thresh_max + pval_thresh_step_size, pval_thresh_step_size)
        print(p_thresh_list)
        S_o = 0 # start index in sorted list
        S_r = 0 # start index in sorted list
        for p_thresh in p_thresh_list:
            # S_o = sum(p_thresh < pv_nnl for pv_nnl in pvals_observed) # slow
            while S_o < len(pvals_observed) and pvals_observed[S_o] > p_thresh:
                S_o += 1

            # Get number of significant interactions for randomized P-values
            # S_r = sum(p_thresh < pv_nnl for pv_nnl in pvals_randomized) # slow
            while S_r < len(pvals_randomized) and pvals_randomized[S_r] > p_thresh:
                S_r += 1

            # Estimate FDR
            fdr = S_r / S_o
            print("X\t" + "\t" + str(fdr) + "\t" + str(p_thresh) + "\t" + str(S_r) + "\t" + str(S_o))

    def perform_randomization_analysis(self, nominal_alpha: float = 0.01, n_iter: int = 1000):

        # Check input parameters
        if nominal_alpha <= 0 or 1.0 < nominal_alpha:
            raise ValueError("Nominal alpha must be in ]0,1]!")

        # Determine number of significant interactions at nominal alpha

        # Get dictionary that stores the numbers of interactions with n read pairs

        # Use dictionary to get list of P-values for randomized data and determine number of significant interactions

    def _get_rp_inter_dict(self):

        rp_inter_dict = {}
        for d_inter in self._interaction_set.interaction_list:
            if d_inter.rp_total in rp_inter_dict:
                rp_inter_dict[d_inter.rp_total] += 1
            else:
                rp_inter_dict[d_inter.rp_total] = 1
        return rp_inter_dict


    def _get_list_of_p_values_from_randomized_data(self, rp_inter_dict: dict = None):

        random_pval_list = []

        for rp_num, i_num in rp_inter_dict.items():

            # Generate random simple read pair counts for current read pair number
            simple_count_list = list(binom.rvs(rp_num, p=0.5, size=i_num))

            # Calculate P-values and append to list
            for simple_count in simple_count_list:
                pval = self._interaction_set._p_values.get_binomial_nnl_p_value(simple_count, rp_num-simple_count)
                random_pval_list.append(pval)

        return random_pval_list
