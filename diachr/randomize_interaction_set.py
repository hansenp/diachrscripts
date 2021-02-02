from scipy.stats import binom
from numpy import arange, log, random, mean, std
from math import floor, ceil
import matplotlib.pyplot as plt
import warnings
import multiprocessing as mp
import itertools
import time
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

warnings.filterwarnings('always')


class RandomizeInteractionSet:

    def __init__(self, random_seed: int = None, interaction_set: DiachromaticInteractionSet = None):

        # Set random seed to be able to reproduce results
        self._random_seed = random_seed
        if random_seed is None:
            self._random_seed = int(time.time())

        # Interaction set on which analyzes are performed
        self._interaction_set = interaction_set

        # Dictionary with all information on the last FDR procedure performed
        self._fdr_info_dict = None

        # Dictionary with all information on the last randomization performed
        self._randomization_info_dict = None

        # Issue a warning, if the FDR procedure was performed on an interaction set with fewer interactions
        self._FDR_WARN_INPUT_INTER_NUM = 16000

        # Issue a warning, if the number of interactions with enough read pairs required for significance at the
        # determined P-value threshold is smaller
        self._FDR_WARN_MIN_RP_INTER_NUM = 10000

        # Issue a warning, if the number of randomized significant interactions at the determined P-value threshold is
        # smaller
        self._FDR_WARN_RAND_SIG_INTER_NUM = 20

    def get_pval_thresh_at_chosen_fdr_thresh(self,
                                             chosen_fdr_thresh: float = 0.05,
                                             pval_thresh_max: float = 0.05,
                                             pval_thresh_step_size: float = 0.00025,
                                             verbose: bool = False
                                             ):
        """
        This function implements the FDR procedure.

        :param chosen_fdr_thresh: Predefined FDR threshold
        :param pval_thresh_max: Largest P-value threshold
        :param pval_thresh_step_size: Step size to increase the P-value threshold
        :param verbose: If true, messages about progress will be written to the screen
        :return: If true, messages about progress will be written to the screen
        """

        if verbose:
            print("[INFO] Performing FDR procedure ...")

        # Check input parameters
        if chosen_fdr_thresh <= 0 or 1.0 < chosen_fdr_thresh:
            raise ValueError("FDR threshold must be in ]0,1]!")
        if pval_thresh_max <= 0 or 1.0 < pval_thresh_max:
            raise ValueError("Maximum P-value threshold must be in ]0,1]!")
        if pval_thresh_max <= pval_thresh_step_size:
            raise ValueError("Step size must be smaller than maximum P-value threshold!")

        if verbose:
            print("\t[INFO] Getting list of observed P-values ...")

        pvals_observed = []
        for d_inter in self._interaction_set.interaction_list:
            nnl_pval = self._interaction_set._p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)
            pvals_observed.append(nnl_pval)

        if verbose:
            print("\t[INFO] Getting list of randomized P-values ...")

        # Get dictionary that stores the numbers of interactions with n read pairs and list of random P-values
        rp_inter_dict = self._get_rp_inter_dict()
        pvals_randomized = self._get_list_of_p_values_from_randomized_data(rp_inter_dict=rp_inter_dict,
                                                                           random_seed=self._random_seed)

        if verbose:
            print("\t[INFO] Going through list of P-value thresholds and estimate FDR ...")

        # Sort lists of P-values
        pvals_observed.sort(key=float, reverse=True)
        pvals_randomized.sort(key=float, reverse=True)

        # Create list of ascending P-value thresholds and prepare lists for results
        pval_thresh_list = arange(pval_thresh_step_size, pval_thresh_max + pval_thresh_step_size, pval_thresh_step_size)
        nnl_pval_thresh_list = []
        min_rp_num_list = []
        min_rp_num_pval_list = []
        sig_num_r_list = []
        sig_num_o_list = []
        fdr_list = []

        # Set start indices in sorted lists (index corresponds to number of significant interactions)
        sig_num_o = 0
        sig_num_r = 0

        # Go through list of ascending P-value thresholds
        for p_thresh in pval_thresh_list:

            # Get minimum number of read pairs required for significance at current  threshold
            min_rp_num, min_rp_num_pval = self._interaction_set._p_values.find_smallest_significant_n(p_thresh)
            min_rp_num_list.append(min_rp_num)
            min_rp_num_pval_list.append(min_rp_num_pval)

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

        if verbose:
            print("\t[INFO] Looking for largest P-value threshold for which the FDR is below the chosen threshold ...")

        # Find the largest P-value threshold for which the FDR is below the chosen threshold
        idx = 0
        while idx < len(fdr_list) and fdr_list[idx] < chosen_fdr_thresh:
            idx += 1

        # Look for irregularities and issue warnings if necessary
        issued_warnings = [0,0,0]
        if len(self._interaction_set.interaction_list) < self._FDR_WARN_INPUT_INTER_NUM:
            issued_warnings[0] = 1
            warnings.warn('Only ' + str(len(self._interaction_set.interaction_list)) + ' input interactions! ' +
                          'Probably not enough to estimate the FDR.')

        min_rp_num_inter_num = self._interaction_set.get_num_of_inter_with_as_many_or_more_read_pairs(
            min_rp_num_list[idx - 1])
        if min_rp_num_inter_num < self._FDR_WARN_MIN_RP_INTER_NUM:
            issued_warnings[1] = 1
            warnings.warn('Only ' + str(min_rp_num_inter_num) +
                          ' interactions with minimum number of ' + str(min_rp_num_list[idx - 1]) +
                          ' read pairs required for significance at the determined threshold! ' +
                          'Probably not enough to estimate the FDR.')

        if sig_num_r_list[idx - 1] < self._FDR_WARN_RAND_SIG_INTER_NUM:
            issued_warnings[2] = 1
            warnings.warn('Only ' + str(sig_num_r_list[idx - 1]) + ' significant interactions after randomization ' +
                          'at determined threshold. ' +
                          'Probably not enough to estimate the FDR.')

        if verbose:
            print("[INFO] ... done.")

        # Prepare and return dictionary for report
        report_dict = {
            'INPUT_PARAMETERS':
                {
                    'CHOSEN_FDR_THRESH': [chosen_fdr_thresh],
                    'PVAL_THRESH_MAX': [pval_thresh_max],
                    'PVAL_THRESH_STEP_SIZE': [pval_thresh_step_size],
                    'INPUT_INTERACTIONS_NUM': [len(self._interaction_set.interaction_list)],
                    'RANDOM_SEED': [self._random_seed]
                },
            'RESULTS_TABLE':
                {
                    'PVAL_THRESH': pval_thresh_list,
                    'NNL_PVAL_THRESH': nnl_pval_thresh_list,
                    'MIN_RP_NUM': min_rp_num_list,
                    'MIN_RP_NUM_PVAL': min_rp_num_pval_list,
                    'SIG_NUM_R': sig_num_r_list,
                    'SIG_NUM_O': sig_num_o_list,
                    'FDR': fdr_list
                },
            'RESULT_INDEX': [idx - 1],
            'MIN_RP_NUM_INTER_NUM': min_rp_num_inter_num,
            'WARNINGS': issued_warnings
        }
        self._fdr_info_dict = report_dict
        return report_dict

    def get_fdr_info_report(self):
        """
        :return: String that contains information about the last FDR procedure performed.
        """

        report = "[INFO] Report on FDR procedure:" + '\n'
        report += "\t[INFO] Input parameters:" + '\n'
        report += "\t\t[INFO] Chosen FDR threshold: " \
                  + "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['CHOSEN_FDR_THRESH'][0]) + '\n'
        report += "\t\t[INFO] Maximum P-value threshold: " \
                  + "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_MAX'][0]) + '\n'
        report += "\t\t[INFO] P-value threshold step size: " \
                  + "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_STEP_SIZE'][0]) + '\n'
        report += "\t\t[INFO] Total number of interactions: " \
                  + "{:,}".format(self._fdr_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]) + '\n'
        report += "\t\t[INFO] Random seed: " \
                  + str(self._fdr_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]) + '\n'

        result_index = self._fdr_info_dict['RESULT_INDEX'][0]

        report += "\t[INFO] Results:" + '\n'
        report += "\t\t[INFO] Determined P-value threshold: " \
                  + "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH'][result_index]) + '\n'
        report += "\t\t[INFO] Determined -ln(P-value threshold): " \
                      + "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['NNL_PVAL_THRESH'][result_index]) + '\n'
        report += "\t\t[INFO] Minimum read pair number: " \
                  + str(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]) + '\n'
        report += "\t\t[INFO] Smallest possible P-value with " \
                  + str(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]) + " read pairs: " \
                  + "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL'][result_index]) + '\n'
        report += "\t\t[INFO] Number of interactions with " \
                  + str(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]) + " or more read pairs: " \
                      + "{:,}".format(self._fdr_info_dict['MIN_RP_NUM_INTER_NUM']) + '\n'
        report += "\t\t[INFO] Number of significant interactions: " \
                  + "{:,}".format(self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O'][result_index]) + '\n'
        report += "\t\t[INFO] Number of randomized significant interactions: " \
                  + "{:,}".format(self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]) + '\n'
        report += "\t\t[INFO] Estimated FDR: " \
                  + "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['FDR'][result_index]) + '\n'

        if 0 < sum(self._fdr_info_dict['WARNINGS']):
            report += "\t[INFO] Warnings: " + '\n'

            if self._fdr_info_dict['WARNINGS'][0] == 1:
                report += "\t\t[INFO] Not enough input interactions!" + '\n'
                report += "\t\t\t[INFO] Only "\
                        + "{:,}".format(self._fdr_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0])\
                        + ' of ' + "{:,}".format(self._FDR_WARN_INPUT_INTER_NUM) +  '\n'

            if self._fdr_info_dict['WARNINGS'][1] == 1:
                report += "\t\t[INFO] Not enough interactions with enough read pairs required for significance!"  + '\n'
                report += "\t\t\t[INFO] Only " \
                          + "{:,}".format(self._fdr_info_dict['MIN_RP_NUM_INTER_NUM']) \
                          + ' of ' + "{:,}".format(self._FDR_WARN_MIN_RP_INTER_NUM) + '\n'

            if self._fdr_info_dict['WARNINGS'][2] == 1:
                report += "\t\t[INFO] Not enough randomized significant interactions!"  + '\n'
                report += "\t\t\t[INFO] Only " \
                          + "{:,}".format(self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]) \
                          + ' of ' + "{:,}".format(self._FDR_WARN_RAND_SIG_INTER_NUM) + '\n'

        report += "[INFO] End of report." + '\n'

        return report

    def get_fdr_info_table_row(self, out_prefix: str = None):
        """
        :return: String consisting of a header line and a line with values relating to last performed FDR procedure
        """

        # Table tag and prefix for output
        table_row = ":TR_FDR:" + '\t'
        table_row += "OUT_PREFIX" + '\t'

        # Input parameters for FDR procedure
        table_row += "CHOSEN_FDR_THRESH" + '\t'
        table_row += "PVAL_THRESH_MAX" + '\t'
        table_row += "PVAL_THRESH_STEP_SIZE" + '\t'
        table_row += "INPUT_INTER_NUM" + '\t'
        table_row += "RANDOM_SEED" + '\t'

        # Results
        table_row += "PVAL_THRESH" + '\t'
        table_row += "NNL_PVAL_THRESH" + '\t'
        table_row += "MIN_RP_NUM" + '\t'
        table_row += "MIN_RP_NUM_PVAL" + '\t'
        table_row += "MIN_RP_NUM_INTER_NUM" + '\t'
        table_row += "SIG_NUM_R" + '\t'
        table_row += "SIG_NUM_O" + '\t'
        table_row += "FDR" + '\t'

        # Warnings
        table_row += "WARNINGS" + '\n'


        # Table tag and prefix for output
        table_row += ":TR_FDR:" + '\t'
        table_row += out_prefix + '\t'

        # Input parameters for FDR procedure
        table_row += "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['CHOSEN_FDR_THRESH'][0]) + '\t'
        table_row += "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_MAX'][0]) + '\t'
        table_row += "{:.5f}".format(self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_STEP_SIZE'][0]) + '\t'
        table_row += str(self._fdr_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]) + '\t'
        table_row += str(self._fdr_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]) + '\t'

        # Results
        result_index = self._fdr_info_dict['RESULT_INDEX'][0]
        table_row += "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH'][result_index]) + '\t'
        table_row += "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['NNL_PVAL_THRESH'][result_index]) + '\t'
        table_row += str(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]) + '\t'
        table_row += "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL'][result_index]) + '\t'
        table_row += str(self._fdr_info_dict['MIN_RP_NUM_INTER_NUM']) + '\t'
        table_row += str(self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]) + '\t'
        table_row += str(self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O'][result_index]) + '\t'
        table_row += "{:.5f}".format(self._fdr_info_dict['RESULTS_TABLE']['FDR'][result_index]) + '\t'

        # Warnings
        table_row += str(self._fdr_info_dict['WARNINGS'][0]) + ','
        table_row += str(self._fdr_info_dict['WARNINGS'][1]) + ','
        table_row += str(self._fdr_info_dict['WARNINGS'][2]) + '\n'

        return table_row

    def perform_randomization_analysis(self, nominal_alpha: float = 0.01, iter_num: int = 1000, thread_num: int = 0,
                                       verbose: bool = False):
        """
        This function implements the entire randomization analysis.

        :param nominal_alpha: Interactions with a larger P-value are classified as significant
        :param iter_num: Number of iterations performed
        :param thread_num: Number of parallel processes in which the iterations are preformed
        :param verbose: If true, messages about progress will be written to the screen
        """

        if verbose:
            print("[INFO] Performing randomization analysis with " + str(iter_num) + " iterations ...")

        # Check input parameters
        if nominal_alpha <= 0 or 1.0 < nominal_alpha:
            raise ValueError("Nominal alpha must be in ]0,1]!")

        # Get negative natural logarithm of nominal alpha
        nnl_nominal_alpha = -log(nominal_alpha)

        if verbose:
            print("\t[INFO] Determining number of significant interactions at nominal alpha ...")

        sig_num_o = 0
        for d_inter in self._interaction_set.interaction_list:
            nnl_pval = self._interaction_set._p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)
            if nnl_nominal_alpha < nnl_pval:
                sig_num_o += 1

        if verbose:
            print("\t[INFO] Randomizing interactions ...")

        # Get dictionary that stores the numbers of interactions with n read pairs
        min_rp_num, max_p_val = self._interaction_set._p_values.find_smallest_significant_n(nominal_alpha)
        rp_inter_dict = self._get_rp_inter_dict(min_rp_num=min_rp_num)
        i_num_randomized = 0 # Determine number of randomized interactions
        for i_num in rp_inter_dict.values():
            i_num_randomized += i_num

        # Use dictionary to get list of P-values for randomized data and determine number of significant interactions
        args_dict = {
            'RP_INTER_DICT': rp_inter_dict,
            'NNL_PVAL_THRESH': nnl_nominal_alpha,
            'ITER_IDX_RANGE': [],
            'VERBOSE': verbose
        }

        # Perform randomization without or with multiprocessing package
        if thread_num == 0:
            args_dict['ITER_IDX_RANGE'] = range(0, iter_num)
            sig_num_r_list = self._perform_n_iterations(args_dict)
        else:

            # Perform permutation for 'thread_num' batches with balanced numbers of iterations
            pool = mp.Pool(processes=thread_num)

            # Divide iterations into 'thread_num' batches of equal size
            batch_iter_nums = list(itertools.repeat(int(iter_num / thread_num), thread_num))
            if 0 < iter_num - thread_num * int(iter_num / thread_num):
                # If there is a remainder, add it to the first element of list
                batch_iter_nums[0] = batch_iter_nums[0] + iter_num - thread_num * int(iter_num / thread_num)

            # Process batches in a separate process
            results = []
            iter_start_idx = 0
            for batch_iter_num in batch_iter_nums:
                args_dict['ITER_IDX_RANGE'] = range(iter_start_idx, iter_start_idx + batch_iter_num)
                result = pool.apply_async(self._perform_n_iterations, args=(args_dict,))
                results.append(result)
                time.sleep(1)
                iter_start_idx += batch_iter_num

            # Wait until all processes are finished
            [result.wait() for result in results]

            # Combine results from different processes
            batch_results = [p.get() for p in results]
            sig_num_r_list = list(itertools.chain.from_iterable(batch_results))

            # Shut down the pool
            pool.close()

        if verbose:
            print("\t[INFO] Calculating summary statistics ...")

        # Determine number of interactions with more significant interactions than originally observed
        sig_num_r_gt_obs = len([sig_num_r for sig_num_r in sig_num_r_list if sig_num_o < sig_num_r])

        # Calculate average number of significant randomized interactions
        sig_num_r_mean = mean(sig_num_r_list)

        # Calculate standard deviation of significant randomized interactions
        sig_num_r_std = std(sig_num_r_list)

        # Calculate Z-score
        z_score = (sig_num_o - sig_num_r_mean) / sig_num_r_std

        # Estimate FDR
        fdr = sig_num_r_mean/sig_num_o

        if verbose:
            print("[INFO] ... done.")

        # Prepare and return dictionary for report
        report_dict = {
            'INPUT_PARAMETERS':
                {
                    'NOMINAL_ALPHA': [nominal_alpha],
                    'ITER_NUM': [iter_num],
                    'INPUT_INTERACTIONS_NUM': [len(self._interaction_set.interaction_list)],
                    'RANDOM_SEED': [self._random_seed]
                },
            'RESULTS':
                {
                    'SIG_NUM_R_LIST': sig_num_r_list,
                    'SUMMARY':
                         {
                            'SIG_NUM_O': [sig_num_o],
                            'I_NUM_RANDOMIZED': [i_num_randomized],
                            'SIG_NUM_R_GT_OBS': [sig_num_r_gt_obs],
                            'SIG_NUM_R_MEAN': [sig_num_r_mean],
                            'SIG_NUM_R_STD': [sig_num_r_std],
                            'Z_SCORE': [z_score],
                            'FDR':[fdr]
                         }
                }
        }
        self._randomization_info_dict = report_dict
        return report_dict

    def get_randomization_info_report(self):
        """
        :return: String that contains information about the last randomization performed.
        """

        report = "[INFO] Report on randomization analysis:" + '\n'

        report += "\t[INFO] Input parameters:" + '\n'
        report += "\t\t[INFO] Number of input interactions: " \
                  + "{:,}".format(self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]) + '\n'
        report += "\t\t[INFO] Nominal alpha: " \
                  + "{:.5f}".format(self._randomization_info_dict['INPUT_PARAMETERS']['NOMINAL_ALPHA'][0]) + '\n'
        report += "\t\t[INFO] Number of iterations: " \
                  + "{:,}".format(self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0]) + '\n'
        report += "\t\t[INFO] Random seed: " \
                  + str(self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]) + '\n'

        report += "\t[INFO] Results:" + '\n'
        report += "\t\t[INFO] Original number of significant interactions: " \
                  + "{:,}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_O'][0]) + '\n'
        report += "\t\t[INFO] Number of randomized interactions: " \
                  + "{:,}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['I_NUM_RANDOMIZED'][0]) + '\n'
        report += "\t\t[INFO] First 10 significant randomized interaction numbers: " + '\n' \
                  + "\t\t\t" + ", ".join(str(i) for i in self._randomization_info_dict['RESULTS']['SIG_NUM_R_LIST'][:10]) \
                  + ", ..." + '\n'
        report += "\t\t[INFO] Iterations with more significant interactions than observed: " \
                  + "{:,}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_GT_OBS'][0]) + '\n'
        report += "\t\t[INFO] Mean number of significant randomized interactions: " \
                  + "{0:,.2f}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]) + '\n'
        report += "\t\t[INFO] Standard deviation of significant randomized interactions: " \
                  + "{:.2f}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_STD'][0]) + '\n'
        report += "\t\t[INFO] Z-score: " \
                  + "{:.2f}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['Z_SCORE'][0]) + '\n'
        report += "\t\t[INFO] Estimated FDR: " \
                  + "{:.5f}".format(self._randomization_info_dict['RESULTS']['SUMMARY']['FDR'][0]) + '\n'

        report += "[INFO] End of report." + '\n'

        return report

    def _get_rp_inter_dict(self, min_rp_num: int = 0):
        """
        Create dictionary required needed for randomization of interactions.

        :param min_rp_num: Interactions with fewer read pairs are not taken into account
        :return: Dictionary that contains the number of interactions (values) for all occurring read pair numbers
        (keys) for the interaction set of this object.
        """

        rp_inter_dict = {}
        for d_inter in self._interaction_set.interaction_list:
            if min_rp_num <= d_inter.rp_total:
                if d_inter.rp_total in rp_inter_dict:
                    rp_inter_dict[d_inter.rp_total] += 1
                else:
                    rp_inter_dict[d_inter.rp_total] = 1
        return rp_inter_dict

    def _get_list_of_p_values_from_randomized_data(self, rp_inter_dict: dict = None, random_seed: int = None):
        """
        This function generates randomized simple and twisted read pair counts for all interactions of the
        interaction set of this object and calculates associated P-values.

        :param rp_inter_dict: Dictionary generated with function '_get_rp_inter_dict'
        :return: List of P-values for randomized interactions
        """

        # Init random generator
        random.seed(random_seed)

        random_pval_list = []

        for rp_num, i_num in rp_inter_dict.items():

            # Generate random simple read pair counts for current read pair number
            simple_count_list = list(binom.rvs(rp_num, p=0.5, size=i_num))

            # Calculate P-values and append to list
            for simple_count in simple_count_list:
                nnl_pval = self._interaction_set._p_values.get_binomial_nnl_p_value(simple_count, rp_num - simple_count)
                random_pval_list.append(nnl_pval)

        return random_pval_list

    def _perform_one_iteration(self, rp_inter_dict: dict = None, nnl_pval_thresh: float = None, random_seed: int = None):
        """
        This function performs a single iteration of the randomization procedure.

        :param rp_inter_dict: Dictionary that contains the number of interactions for different numbers of read pairs
        :param nnl_pval_thresh: Negative natural logarithm of P-value threshold
        :param random_seed: Number used to init random generator
        :return: Number of randomized significant interactions
        """

        # Get list of P-values for randomized interactions
        randomized_nnl_pvals = self._get_list_of_p_values_from_randomized_data(rp_inter_dict, random_seed=random_seed)

        # Determine number of randomized significant interactions at given P-value threshold
        sig_num_r = len([1 for nnl_pval in randomized_nnl_pvals if nnl_pval > nnl_pval_thresh])

        return sig_num_r

    def _perform_n_iterations(self, args_dict: dict = None):
        """
        This function performs a given number of iterations of the randomization procedure.

        :param param_dict: Dictionary containing all function arguments
        :return: List of numbers of randomized significant interactions for each iteration
        """

        # We pass all function parameters in a dictionary, because that's easier to use with 'pool.apply_async'
        rp_inter_dict = args_dict['RP_INTER_DICT']
        nnl_pval_thresh = args_dict['NNL_PVAL_THRESH']
        iter_idx_range = args_dict['ITER_IDX_RANGE']
        verbose = args_dict['VERBOSE']

        # Init list with numbers of randomized significant interactions for each iteration
        sig_num_r_list = []

        if verbose:
            print("\t\t[INFO] Performing " + str(len(iter_idx_range)) + " iterations ...")
            print("\t\t\t[INFO] First iteration indices: " + ", ".join(str(i) for i in iter_idx_range[:10]) + ", ...")

        # Perform each iteration with its own random seed that corresponds by adding the iteration index
        for iter_idx in iter_idx_range:
            sig_num_r_list.append(self._perform_one_iteration(rp_inter_dict=rp_inter_dict,
                                                              nnl_pval_thresh=nnl_pval_thresh,
                                                              random_seed=self._random_seed + iter_idx))

        return sig_num_r_list

    def get_fdr_info_plot(self, pdf_file_name: str = None, analysis_name: str = None):
        """
        This function creates a graphical representation of the results from the last FDR procedure performed.

        :param pdf_file_name: Name of the PDF file that will be created
        :param analysis_name: Name of of the analysis that will be shown in the graphical representation
        :return: Nothing, if no FDR procedure has been performed yet or otherwise a 'Figure' object of 'matplotlib'
        that can be displayed in a Jupyter notebook
        """

        # Check whether an analysis has already been performed
        if self._fdr_info_dict is None:
            print("[ERROR] No analysis has been performed yet! There is nothing to plot.")
            return

        # Set common plot parameters
        hv_lwd = 0.5  # line width of horizontal and vertical red dashed lines
        hv_col = 'red'  # color of horizontal and vertical red dashed lines
        header_font_size = 10

        # Extract data from fdr_info_dict
        # -------------------------------

        # Input parameters
        chosen_fdr_thresh = self._fdr_info_dict['INPUT_PARAMETERS']['CHOSEN_FDR_THRESH'][0]
        pval_thresh_max = self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_MAX'][0]
        pval_thresh_step_size = self._fdr_info_dict['INPUT_PARAMETERS']['PVAL_THRESH_STEP_SIZE'][0]
        input_interactions_num = self._fdr_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]
        random_seed = self._fdr_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]
        result_index = self._fdr_info_dict['RESULT_INDEX'][0]

        # Results
        pval_thresh_column = self._fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH']
        pval_thresh_result = self._fdr_info_dict['RESULTS_TABLE']['PVAL_THRESH'][result_index]
        nnl_pval_thresh_result = self._fdr_info_dict['RESULTS_TABLE']['NNL_PVAL_THRESH'][result_index]
        min_rp_num_column = self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM']
        min_rp_num_result = self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM'][result_index]
        min_rp_num_pval_column = self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL']
        min_rp_num_pval_result = self._fdr_info_dict['RESULTS_TABLE']['MIN_RP_NUM_PVAL'][result_index]
        min_rp_num_inter_num = self._fdr_info_dict['MIN_RP_NUM_INTER_NUM']
        sig_num_o_column = self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O']
        sig_num_o_result = self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_O'][result_index]
        sig_num_r_column = self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R']
        sig_num_r_result = self._fdr_info_dict['RESULTS_TABLE']['SIG_NUM_R'][result_index]
        fdr_column = self._fdr_info_dict['RESULTS_TABLE']['FDR']
        fdr_result = self._fdr_info_dict['RESULTS_TABLE']['FDR'][result_index]
        fdr_max = max(self._fdr_info_dict['RESULTS_TABLE']['FDR'])
        if self._fdr_info_dict['WARNINGS'][0] == 1:
            warn_input_inter_col = 'r'
        else:
            warn_input_inter_col = 'k'
        if self._fdr_info_dict['WARNINGS'][1] == 1:
            warn_min_rp_inter_col = 'r'
        else:
            warn_min_rp_inter_col = 'k'
        if self._fdr_info_dict['WARNINGS'][2] == 1:
            warn_rand_sig_inter_col = 'r'
        else:
            warn_rand_sig_inter_col = 'k'

        # Create figure with plots for individual columns
        # -----------------------------------------------

        fig, ax = plt.subplots(6, figsize=(7, 19), gridspec_kw={'height_ratios': [2, 1, 1, 1, 1, 1]})

        # Add field with information about analysis
        plt.plot(ax=ax[0])
        ax[0].spines['left'].set_color('white')
        ax[0].spines['right'].set_color('white')
        ax[0].spines['top'].set_color('white')
        ax[0].spines['bottom'].set_color('white')
        ax[0].tick_params(axis='x', colors='white')
        ax[0].tick_params(axis='y', colors='white')
        ax[0].text(-0.2, 1.00, 'FDR results', fontsize=18, fontweight='bold')
        ax[0].text(-0.18, 0.90, 'Analysis name: ' + analysis_name, fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.85, 'Chosen FDR threshold: ' + "{:.5f}".format(chosen_fdr_thresh),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.80, 'Maximum P-value threshold: ' + "{:.5f}".format(pval_thresh_max),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.75, 'P-value threshold step size: ' + "{:.5f}".format(pval_thresh_step_size),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.70, 'Total number of interactions: ' + "{:,}".format(input_interactions_num),
                   fontsize=header_font_size, fontweight='bold', color=warn_input_inter_col)
        ax[0].text(-0.18, 0.65, 'Random seed: ' + str(random_seed),
                   fontsize=header_font_size)

        ax[0].text(-0.18, 0.55, 'Determined P-value threshold: ' + "{:.5f}".format(pval_thresh_result),
                   fontsize=header_font_size,
                   fontweight='bold')
        ax[0].text(-0.18, 0.50, 'Determined -ln(P-value threshold): ' + "{:.5f}".format(nnl_pval_thresh_result),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.45, 'Minimum read pair number: ' + str(min_rp_num_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.40,
                   'Smallest possible P-value with ' + str(min_rp_num_result) + ' read pairs: ' + "{:.5f}".format(
                       min_rp_num_pval_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.35,
                   'Number of interactions with ' + str(min_rp_num_result) + ' or more read pairs: ' + "{:,}".format(
                       min_rp_num_inter_num), fontsize=header_font_size, fontweight='bold', color=warn_min_rp_inter_col)
        ax[0].text(-0.18, 0.30, 'Number of significant interactions: ' + "{:,}".format(sig_num_o_result),
                   fontsize=header_font_size,
                   fontweight='bold')
        ax[0].text(-0.18, 0.25, 'Number of randomized significant interactions: ' + "{:,}".format(sig_num_r_result),
                   fontsize=header_font_size, color=warn_rand_sig_inter_col)
        ax[0].text(-0.18, 0.20, 'Estimated FDR: ' + "{:.5f}".format(fdr_result), fontsize=header_font_size,
                   fontweight='bold')

        # Plot P-value thresholds
        ax[1].plot(pval_thresh_column, pval_thresh_column)
        ax[1].set_title('P-value threshold: ' + "{:.5f}".format(pval_thresh_result), loc='left')
        ax[1].set(xlabel='P-value threshold')
        ax[1].set(ylabel='PVAL_TRESH')
        ax[1].axhline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[1].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot minimum read pair numbers
        ax[2].plot(pval_thresh_column, min_rp_num_column)
        ax[2].set_title('Minimum read pair number: ' + str(min_rp_num_result), loc='left')
        ax[2].set(xlabel='P-value threshold')
        ax[2].set(ylabel='MIN_RP_NUM')
        ax[2].axhline(min_rp_num_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[2].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Smallest possible P-values with minimum read pair numbers
        ax[3].plot(pval_thresh_column, min_rp_num_pval_column)
        ax[3].set_title('Smallest possible P-value with ' + str(min_rp_num_result) + ' read pairs: ' + "{:.5f}".format(
            min_rp_num_pval_result), loc='left')
        ax[3].set(xlabel='P-value threshold')
        ax[3].set(ylabel='MIN_RP_NUM_PVAL')
        ax[3].axhline(min_rp_num_pval_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[3].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot number of significant interactions
        ax[4].plot(pval_thresh_column, sig_num_r_column, label='SIG_NUM_R')
        ax[4].plot(pval_thresh_column, sig_num_o_column, label='SIG_NUM_O')
        ax[4].set_title('Number of significant interactions: ' + "{:,}".format(sig_num_o_result) + ' (' + "{:,}".format(
            sig_num_r_result) + ')', loc='left')
        ax[4].set(xlabel='P-value threshold')
        ax[4].legend(fontsize=8)
        ax[4].axhline(sig_num_o_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[4].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot estimated FDR
        ax[5].plot(pval_thresh_column, fdr_column)
        ax[5].set_title('Estimated FDR: ' + "{:.5f}".format(fdr_result), loc='left')
        ax[5].set(xlabel='P-value threshold')
        ax[5].set(ylabel='FDR')
        ax[5].axhline(fdr_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[5].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        ax[5].text(pval_thresh_max - pval_thresh_max/5, fdr_result + fdr_max/16, 'FDR: ' + "{:.5f}".format(fdr_result),
                 bbox={'color': 'lightblue', 'alpha': 0.5, 'pad': 4})

        ax[5].text(pval_thresh_result + pval_thresh_max/60, fdr_result - fdr_max/9, 'P-value: ' + "{:.5f}".format(pval_thresh_result),
                 bbox={'color': 'lightblue', 'alpha': 0.5, 'pad': 4})

        # Finalize save to PDF and return 'Figure' object
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    def get_randomization_info_plot(self, pdf_file_name: str = None, analysis_name: str = None):
        """
        This function creates a graphical representation of the results from the last randomization analysis performed.

        :param pdf_file_name: Name of the PDF file that will be created
        :param analysis_name: Name of of the analysis that will be shown in the graphical representation
        :return: Nothing, if no FDR procedure has been performed yet or otherwise a 'Figure' object of 'matplotlib'
        that can be displayed in a Jupyter notebook
        """

        # Check whether an analysis has already been performed
        if self._randomization_info_dict is None:
            print("[ERROR] No analysis has been performed yet! There is nothing to plot.")
            return

        # Set common plot parameters
        header_font_size = 10
        bin_face_color = 'blue'
        bin_edge_color = 'darkblue'

        # Extract data from fdr_info_dict
        # -------------------------------

        # Input parameters
        pdf_file_name = pdf_file_name
        analysis_name = analysis_name
        nominal_alpha = self._randomization_info_dict['INPUT_PARAMETERS']['NOMINAL_ALPHA'][0]
        iter_num = self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0]
        input_interactions_num = self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]
        random_seed = self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]

        # Results
        sig_num_r_list = self._randomization_info_dict['RESULTS']['SIG_NUM_R_LIST']
        sig_num_o = self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_O'][0]
        i_num_randomized = self._randomization_info_dict['RESULTS']['SUMMARY']['I_NUM_RANDOMIZED'][0]
        sig_num_r_gt_obs = self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_GT_OBS'][0]
        sig_num_r_mean = self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_MEAN'][0]
        sig_num_r_sd = self._randomization_info_dict['RESULTS']['SUMMARY']['SIG_NUM_R_STD'][0]
        z_score = float(self._randomization_info_dict['RESULTS']['SUMMARY']['Z_SCORE'][0])
        fdr = float(self._randomization_info_dict['RESULTS']['SUMMARY']['FDR'][0])

        # Highlight text if there are iterations with more significant interactions than originally observed
        if 0 < sig_num_r_gt_obs:
            sig_num_r_gt_obs_color = 'red'
            sig_num_r_gt_obs_fontweight = 'bold'
        else:
            sig_num_r_gt_obs_color = 'black'
            sig_num_r_gt_obs_fontweight = 'normal'

        # Create figure with text field and two histograms
        # ------------------------------------------------

        fig, ax = plt.subplots(3, figsize=(7, 10.86), gridspec_kw={'height_ratios': [2, 1, 1]})

        # Add field with information about analysis
        plt.plot(ax=ax[0])
        ax[0].spines['left'].set_color('white')
        ax[0].spines['right'].set_color('white')
        ax[0].spines['top'].set_color('white')
        ax[0].spines['bottom'].set_color('white')
        ax[0].tick_params(axis='x', colors='white')
        ax[0].tick_params(axis='y', colors='white')
        ax[0].text(-0.2, 1.00, 'Randomization analysis results', fontsize=18, fontweight='bold')
        ax[0].text(-0.18, 0.90, 'Analysis name: ' + analysis_name, fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.85, 'Number of input interactions: ' + "{:,}".format(input_interactions_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.80, 'Nominal alpha: ' + "{:.5f}".format(nominal_alpha),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.75, 'Number of iterations: ' + "{:,}".format(iter_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.70, 'Random seed: ' + str(random_seed),
                   fontsize=header_font_size)

        ax[0].text(-0.18, 0.60, 'Number of observed significant interactions: ' + "{:,}".format(sig_num_o),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.55, 'Number of randomized interactions: ' + "{:,}".format(i_num_randomized),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.50,
                   'Iterations with more significant interactions than observed: ' + "{:,}".format(sig_num_r_gt_obs),
                   fontsize=header_font_size, color=sig_num_r_gt_obs_color, fontweight=sig_num_r_gt_obs_fontweight)
        ax[0].text(-0.18, 0.45,
                   'Mean number of significant randomized interactions: ' + "{:,.2f}".format(sig_num_r_mean),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.40,
                   'Standard deviation of significant randomized interactions: ' + "{:.2f}".format(sig_num_r_sd),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.35, 'Z-score: ' + "{:.2f}".format(z_score), fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.30, 'Estimated FDR : ' + "{:.5f}".format(fdr), fontsize=header_font_size, fontweight='bold')

        # Plot distribution of iterations centered on mean of randomized significant interaction numbers
        x_lim_left = floor(sig_num_r_mean - 5 * sig_num_r_sd)
        x_lim_right = ceil(sig_num_r_mean + 5 * sig_num_r_sd)
        bin_width = ceil(sig_num_r_sd / 3)
        bins = range(x_lim_left, x_lim_right + bin_width, bin_width)
        n, bins, patches = ax[1].hist(sig_num_r_list, bins=bins, density=False, facecolor=bin_face_color,
                                      edgecolor=bin_edge_color, alpha=0.75)
        ax[1].set_xlim(x_lim_left, x_lim_right)
        ax[1].set_title('Centered on mean of randomized significant interaction numbers', loc='left')
        ax[1].set_xlabel('Randomized significant interactions')
        ax[1].set_ylabel('Iterations')
        ax[1].axvspan(sig_num_r_mean - 3 * sig_num_r_sd, sig_num_r_mean + 3 * sig_num_r_sd, color='lightblue',
                      alpha=0.9, zorder=0)
        ax[1].axvline(sig_num_r_mean, linestyle='--', color='red', linewidth=1)
        ax[1].axvline(sig_num_r_mean - 2 * sig_num_r_sd, linestyle='--', color='w', linewidth=1, zorder=0)
        ax[1].axvline(sig_num_r_mean - 1 * sig_num_r_sd, linestyle='--', color='w', linewidth=1, zorder=0)
        ax[1].axvline(sig_num_r_mean + 1 * sig_num_r_sd, linestyle='--', color='w', linewidth=1, zorder=0)
        ax[1].axvline(sig_num_r_mean + 2 * sig_num_r_sd, linestyle='--', color='w', linewidth=1, zorder=0)
        y_pos_lab = max(n) - max(n) / 15
        x_pos_lab = sig_num_r_mean - 4.75 * sig_num_r_sd
        ax[1].text(x_pos_lab, y_pos_lab, 'Mean: ' + "{:.2f}".format(sig_num_r_mean),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})
        y_pos_lab = max(n) - 3 * (max(n) / 15)
        ax[1].text(x_pos_lab, y_pos_lab, 'Standard deviation: ' + "{:.2f}".format(sig_num_r_sd),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})

        # Plot distribution of iterations in relation to observed number of significant interactions
        x_lim_left = floor(sig_num_r_mean - 5 * sig_num_r_sd)
        x_lim_right = ceil(sig_num_o + 5 * sig_num_r_sd)
        bin_width = ceil(sig_num_r_sd / 3)
        bins = range(x_lim_left, x_lim_right + bin_width, bin_width)
        n, bins, patches = ax[2].hist(sig_num_r_list, bins=bins, density=False, facecolor=bin_face_color,
                                      edgecolor=bin_edge_color, alpha=0.75)
        ax[2].set_xlim(sig_num_r_mean - 5 * sig_num_r_sd, sig_num_o + 5 * sig_num_r_sd)
        ax[2].set_title('In relation to observed number of significant interactions', loc='left')
        ax[2].set_xlabel('Randomized significant interactions')
        ax[2].set_ylabel('Iterations')
        ax[2].axvspan(sig_num_r_mean - 3 * sig_num_r_sd, sig_num_r_mean + 3 * sig_num_r_sd, color='lightblue',
                      alpha=0.9, zorder=0)
        ax[2].axvline(sig_num_r_mean, linestyle='--', color='red', linewidth=1)
        ax[2].axvline(sig_num_o, linestyle='--', color='red', linewidth=1)
        y_pos_lab = max(n) - max(n) / 15
        x_pos_lab = sig_num_r_mean + 6 * sig_num_r_sd
        ax[2].text(x_pos_lab, y_pos_lab, 'Observed: ' + "{:,}".format(sig_num_o),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})
        y_pos_lab = max(n) - 3 * (max(n) / 15)
        ax[2].text(x_pos_lab, y_pos_lab, 'Z-score: ' + "{:.2f}".format(z_score),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})

        # Finalize save to PDF and return 'Figure' object
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig



