from scipy.stats import binom
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

        # Use dictionary to get list of P-values for randomized data

        # Determine numbers of significant interactions in observed and randomized data for increasing P-value
        # thresholds and estimate FDR.

    def perform_randomization_analysis(self, nominal_alpha: float = 0.01, n_iter: int = 1000):

        # Check input parameters
        if nominal_alpha <= 0 or 1.0 < nominal_alpha:
            raise ValueError("Nominal alpha must be in ]0,1]!")

        # Determine number of significant interactions at nominal alpha

        # Get dictionary that stores the numbers of interactions with n read pairs

        # Use dictionary to get list of P-values for randomized data and determine number of significant interactions

    def _get_rp_inter_dict(self):

        rp_inter_dict = {}


    def _get_list_of_p_values_from_randomized_data(self, rp_inter_dict: dict = None):

        pval_list = []

        for n, i_num in self._N_DICT.items():

            # Generate random simple read pair counts for current n
            simple_count_list = list(binom.rvs(n, p=0.5, size=i_num))

            # Calculate P-values and append to list
            for simple_count in simple_count_list:
                pval = self._interaction_set._p_values.get_binomial_nnl_p_value(simple_count, n-simple_count)
                pval_list.append(pval)

        return pval_list
