from scipy.stats import binom
from numpy import arange, log, random
import warnings
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
warnings.filterwarnings('always')

class RandomizeInteractionSet:

    def __init__(self, random_seed: int = 42, interaction_set: DiachromaticInteractionSet = None):
        self._random_seed = random_seed
        random.seed(random_seed)
        self._interaction_set = interaction_set

    def get_pval_tresh_at_chosen_fdr_tresh(self,
                                           chosen_fdr_thresh: float = 0.05,
                                           pval_thresh_max: float = 0.05,
                                           pval_thresh_step_size: float = 0.00025
                                           ):

        # Set random seed
        random.seed(self._random_seed)

        # Check input parameters
        if chosen_fdr_thresh <= 0 or 1.0 < chosen_fdr_thresh:
            raise ValueError("FDR threshold must be in ]0,1]!")
        if pval_thresh_max <= 0 or 1.0 < pval_thresh_max:
            raise ValueError("Maximum P-value threshold must be in ]0,1]!")
        if pval_thresh_max <= pval_thresh_step_size:
            raise ValueError("Step size must be smaller than maximum P-value threshold!")

        # Get list of observed P-values
        pvals_observed = []
        for d_inter in self._interaction_set.interaction_list:
            pvals_observed.append(d_inter._nln_pval)

        # Get dictionary that stores the numbers of interactions with n read pairs and list of random P-values
        rp_inter_dict = self._get_rp_inter_dict()
        pvals_randomized = self._get_list_of_p_values_from_randomized_data(rp_inter_dict=rp_inter_dict)

        # Sort lists of P-values
        pvals_observed.sort(key=float, reverse=True)
        pvals_randomized.sort(key=float, reverse=True)

        # Create list of ascending P-value thresholds and prepare lists for results
        pval_thresh_list = arange(pval_thresh_step_size, pval_thresh_max + pval_thresh_step_size, pval_thresh_step_size)
        nnl_pval_thresh_list = []
        sig_num_r_list = []
        sig_num_o_list = []
        fdr_list = []

        # Set start indices in sorted lists (index corresponds to number of significant interactions)
        sig_num_o = 0
        sig_num_r = 0

        # Go through list of ascending P-value thresholds
        for p_thresh in pval_thresh_list:

            # Get negative natural logarithm of P-value threshold
            nnl_pval_thresh = -log(p_thresh)
            nnl_pval_thresh_list.append(nnl_pval_thresh)

            # Determine number of observed significant interactions for current threshold
            while sig_num_o < len(pvals_observed) and pvals_observed[sig_num_o] > nnl_pval_thresh:
                sig_num_o += 1
            sig_num_o_list.append(sig_num_o)

            # Determine number of randomized significant interactions for current threshold
            while sig_num_r < len(pvals_randomized) and pvals_randomized[sig_num_r] > nnl_pval_thresh:
                sig_num_r += 1
            sig_num_r_list.append(sig_num_r)

            # Estimate FDR
            fdr = sig_num_r / sig_num_o
            fdr_list.append(fdr)

        # Check whether FDR increases monotonically
        if all(x<=y for x,y in zip(fdr_list,fdr_list[1:])):
            warnings.warn("FDR does not grow monotonically with P-value threshold!" + '\n' + \
                          "Take a close look at the results for each of the P-value thresholds.")

        # Find the largest P-value threshold for which the FDR is below the chosen threshold
        idx = 0
        while idx < len(fdr_list) and fdr_list[idx] < chosen_fdr_thresh:
            idx += 1

        # Prepare and return dictionary for report
        report_dict = {
            'RESULTS_TABLE':
                {
                'P_VAL_TRESH': pval_thresh_list,
                'NNL_P_VAL_TRESH': nnl_pval_thresh_list,
                'SIG_NUM_R': sig_num_r_list,
                'SIG_NUM_O': sig_num_o_list,
                'FDR': fdr_list
                },
            'ROW_INDEX': [idx-1]
        }
        return report_dict

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
