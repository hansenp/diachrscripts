from scipy.stats import binom
from numpy import arange, log, random, mean, std
from math import floor, ceil
import matplotlib.pyplot as plt
import warnings
import multiprocessing as mp
import itertools
import time
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.diachromatic_interaction_set import BinomialModel

warnings.filterwarnings('always')


class RandomizeInteractionSet:

    def __init__(self, random_seed: int = None):

        # Set random seed to be able to reproduce results
        self._random_seed = random_seed
        if random_seed is None:
            self._random_seed = int(time.time())

        # Read pair interaction dictionary used to get random simple read pair counts
        self._rp_inter_dict = None

        # Object for the calculation and storage of already calculated P-values
        self.p_values = BinomialModel()

        # Negative natural logarithm of nominal alpha
        self._nnl_nominal_alpha = None

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
                                             interaction_set: DiachromaticInteractionSet = None,
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
        for d_inter in interaction_set.interaction_list:
            nnl_pval = self.p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)
            pvals_observed.append(nnl_pval)

        if verbose:
            print("\t[INFO] Getting list of randomized P-values ...")

        # Get dictionary that stores the numbers of interactions with n read pairs and list of random P-values
        self._rp_inter_dict = self._get_rp_inter_dict(interaction_set=interaction_set)
        pvals_randomized = self._get_list_of_p_values_from_randomized_data(interaction_set=interaction_set,
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
            min_rp_num, min_rp_num_pval = interaction_set._p_values.find_smallest_significant_n(p_thresh)
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
        issued_warnings = [0, 0, 0]
        if len(interaction_set.interaction_list) < self._FDR_WARN_INPUT_INTER_NUM:
            issued_warnings[0] = 1
            warnings.warn('Only ' + str(len(interaction_set.interaction_list)) + ' input interactions! ' +
                          'Probably not enough to estimate the FDR.')

        min_rp_num_inter_num = interaction_set.get_num_of_inter_with_as_many_or_more_read_pairs(
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
                    'INPUT_INTERACTIONS_NUM': [len(interaction_set.interaction_list)],
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
                report += "\t\t\t[INFO] Only " \
                          + "{:,}".format(self._fdr_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]) \
                          + ' of ' + "{:,}".format(self._FDR_WARN_INPUT_INTER_NUM) + '\n'

            if self._fdr_info_dict['WARNINGS'][1] == 1:
                report += "\t\t[INFO] Not enough interactions with enough read pairs required for significance!" + '\n'
                report += "\t\t\t[INFO] Only " \
                          + "{:,}".format(self._fdr_info_dict['MIN_RP_NUM_INTER_NUM']) \
                          + ' of ' + "{:,}".format(self._FDR_WARN_MIN_RP_INTER_NUM) + '\n'

            if self._fdr_info_dict['WARNINGS'][2] == 1:
                report += "\t\t[INFO] Not enough randomized significant interactions!" + '\n'
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

    def perform_randomization_analysis(self,
                                       interaction_set: DiachromaticInteractionSet = None,
                                       nominal_alphas: [float] = [0.01],
                                       iter_num: int = 1000,
                                       thread_num: int = 0,
                                       verbose: bool = False):
        """
        This function implements the entire randomization analysis.

        :param interaction_set: Input interaction set
        :param nominal_alphas: List of nominal alphas
        :param iter_num: Number of iterations that will be performed
        :param thread_num: Number of parallel processes in which batches iterations are preformed
        :param verbose: If true, messages about progress will be written to the screen
        """

        if verbose:
            print("[INFO] Performing randomization analysis with " + str(iter_num) + " iterations ...")

        # Check input parameters
        for nominal_alpha in nominal_alphas:
            if nominal_alpha <= 0 or 1.0 < nominal_alpha:
                raise ValueError("Nominal alpha must be in ]0,1]!")
        nominal_alphas = sorted(nominal_alphas)
        nominal_alphas = [float("{:.5f}".format(i)) for i in nominal_alphas]

        # Get negative natural logarithm of nominal alpha
        nnl_nominal_alphas = -log(nominal_alphas)

        if verbose:
            print("\t[INFO] Determining number of significant interactions at each nominal alpha ...")

        # Get list of observed P-values
        nnl_pval_list = []
        for d_inter in interaction_set.interaction_list:
            nnl_pval_list.append(self.p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted))

        # Determine number of significant interactions for each nominal alpha
        sig_num_o_list = self._determine_significant_pvals_at_nominal_alphas_nnl(nnl_nominal_alphas=nnl_nominal_alphas, nnl_p_values=nnl_pval_list)

        if verbose:
            print("\t[INFO] Randomizing interactions ...")

        # Get dictionary that stores the numbers of interactions with n read pairs
        min_rp_num, max_p_val = interaction_set._p_values.find_smallest_significant_n(max(nominal_alphas))
        self._rp_inter_dict = self._get_rp_inter_dict(interaction_set=interaction_set, min_rp_num=min_rp_num)
        i_num_randomized = 0  # Determine number of randomized interactions
        for i_num in self._rp_inter_dict.values():
            i_num_randomized += i_num

        if verbose:
            print("\t\t[INFO] Created RP_INTER_DICT for " + "{:,}".format(i_num_randomized) + " interactions ...")

        # Perform randomization without or with multiprocessing package
        if thread_num == 0:
            iter_start_idx = 0
            sig_num_r_list_of_lists = self._perform_n_iterations(iter_start_idx, iter_num, nnl_nominal_alphas, verbose)
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
                result = pool.apply_async(self._perform_n_iterations, args=(iter_start_idx, iter_start_idx + batch_iter_num, nnl_nominal_alphas, verbose))
                results.append(result)
                iter_start_idx += batch_iter_num

            # Wait until all processes are finished
            [result.wait() for result in results]

            # Combine results from different processes
            batch_results = [p.get() for p in results]
            sig_num_r_list_of_lists = list(itertools.chain.from_iterable(batch_results))

            # Shut down the pool
            pool.close()

        if verbose:
            print("\t[INFO] Combining results from all iterations for different nominal alphas ...")

        sig_num_r_lists_list = {}
        sig_num_r_mean_list = []
        sig_num_r_sd_list = []
        sig_num_r_gt_obs_list = []
        z_score_list = []
        fdr_list = []
        issue_z_score_warning = False
        for nominal_alpha_idx in range(0, len(nominal_alphas)):

            # Get current nominal alpha
            nominal_alpha = nominal_alphas[nominal_alpha_idx]
            sig_num_r_lists_list[nominal_alpha] = []

            # Extract numbers of significant randomized interactions for current nominal alpha
            for iter_idx in range(0, len(sig_num_r_list_of_lists)):
                sig_num_r_lists_list[nominal_alpha].append(sig_num_r_list_of_lists[iter_idx][nominal_alpha_idx])

            # Collect and calculate results for current alpha
            sig_num_r_mean = mean(sig_num_r_lists_list[nominal_alpha])
            sig_num_r_sd = std(sig_num_r_lists_list[nominal_alpha])
            sig_num_o = sig_num_o_list[nominal_alpha_idx]
            sig_num_r_gt_obs = len([sig_num_r for sig_num_r in sig_num_r_lists_list[nominal_alpha] if sig_num_o < sig_num_r])
            if sig_num_r_sd != 0:
                z_score = (sig_num_o - sig_num_r_mean) / sig_num_r_sd
            else:
                z_score = "NA"
                issue_z_score_warning = True
            fdr = sig_num_r_mean / sig_num_o

            # Append results to list
            sig_num_r_mean_list.append(sig_num_r_mean)
            sig_num_r_sd_list.append(sig_num_r_sd)
            sig_num_r_gt_obs_list.append(sig_num_r_gt_obs)
            z_score_list.append(z_score)
            fdr_list.append(fdr)

        if verbose:
            print("[INFO] ... done.")

        # Report warnings
        if issue_z_score_warning:
            warnings.warn('Failed to calculate a Z-score for at least one nominal alpha!')

        # Prepare and return dictionary for report
        report_dict = {
            'INPUT_PARAMETERS':
                {
                    'ITER_NUM': [iter_num],
                    'INPUT_INTERACTIONS_NUM': [len(interaction_set.interaction_list)],
                    'RANDOM_SEED': [self._random_seed],
                },
            'SIG_NUM_R_LISTS': sig_num_r_lists_list,
            'I_NUM_RANDOMIZED': i_num_randomized,
            'RESULTS':
                {
                    'NOMINAL_ALPHA': nominal_alphas,
                    'SIG_NUM_R_GT_OBS': sig_num_r_gt_obs_list,
                    'SIG_NUM_O': sig_num_o_list,
                    'SIG_NUM_R_MEAN': sig_num_r_mean_list,
                    'SIG_NUM_R_SD': sig_num_r_sd_list,
                    'Z_SCORE': z_score_list,
                    'FDR': fdr_list
                }
        }
        self._randomization_info_dict = report_dict
        return report_dict

    def get_largest_nominal_alpha_index_at_chosen_fdr_thresh(self, chosen_fdr_thresh: float = 0.05):
        """
        This function determines the index of the largest nominal alpha for which the FDR is still below the passed
        threshold.

        :param chosen_fdr_thresh: FDR threshold
        :return: Index of the largest nominal alpha for which the FDR is still below the passed threshold
        """

        # Check whether an analysis has already been performed
        if self._randomization_info_dict is None:
            print("[ERROR] No analysis has been performed yet! Cannot determine P-value threshold.")
            return

        # Determine results index for chosen FDR threshold
        result_index = 0
        for fdr_idx in range(0, len(self._randomization_info_dict['RESULTS']['FDR'])):
            if self._randomization_info_dict['RESULTS']['FDR'][fdr_idx] <= chosen_fdr_thresh:
                result_index = fdr_idx

        return result_index

    def get_randomization_info_report_at_chosen_fdr_threshold(self, chosen_fdr_threshold: float = 0.05):

        result_index = self.get_largest_nominal_alpha_index_at_chosen_fdr_thresh(chosen_fdr_threshold)
        p_thresh = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA'][result_index]
        info_report = self.get_randomization_info_report(p_thresh)

        return info_report

    def get_randomization_info_report(self, nominal_alpha: float = None):
        """
        :return: String that contains information about the last randomization performed.
        """

        # Check if there are any results for the passed nominal alpha
        nominal_alphas = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']
        if nominal_alpha not in nominal_alphas:
            if nominal_alpha is None:
                report = "[ERROR] No nominal alpha was passed!" + '\n'
            else:
                report = "[ERROR] No results available for a nominal alpha of " + "{:.5f}".format(nominal_alpha) + '!' + '\n'
            report += "\t[ERROR] Results available for: " + ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas) +  '\n'
            return report
        nominal_alpha_idx = nominal_alphas.index(nominal_alpha)

        report = "[INFO] Report on randomization:" + '\n'

        report += "\t[INFO] Input parameters:" + '\n'
        report += "\t\t[INFO] Number of input interactions: " \
                  + "{:,}".format(self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]) + '\n'
        report += "\t\t[INFO] Number of iterations: " \
                  + "{:,}".format(self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0]) + '\n'
        report += "\t\t[INFO] Random seed: " \
                  + str(self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]) + '\n'
        report += "\t\t[INFO] Available nominal alphas: " + '\n'
        if len(nominal_alphas) < 5:
            report += '\t\t\t' + ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas) + '\n'
        else:
            report += '\t\t\t' + ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas[:3]) + ", ..., " + "{:.5f}".format(nominal_alphas[-1]) +'\n'

        report += "\t[INFO] Results for a nominal alpha of " + "{:.5f}".format(nominal_alphas[nominal_alpha_idx]) + ":" + '\n'
        report += "\t\t[INFO] Number of randomized interactions: " \
                  + "{:,}".format(self._randomization_info_dict['I_NUM_RANDOMIZED']) + '\n'
        report += "\t\t[INFO] Significant randomized interaction numbers: " + '\n'
        if len(self._randomization_info_dict['SIG_NUM_R_LISTS'][nominal_alphas[nominal_alpha_idx]]) < 10:
            report += '\t\t\t' + ", ".join(str(i) for i in self._randomization_info_dict['SIG_NUM_R_LISTS'][nominal_alphas[nominal_alpha_idx]]) +  '\n'
        else:
            report += '\t\t\t' + ", ".join(str(i) for i in self._randomization_info_dict['SIG_NUM_R_LISTS'][nominal_alphas[nominal_alpha_idx]][:9]) + ", ..." + '\n'
        report += "\t\t[INFO] Iterations with more significant interactions than observed: " \
                  + "{:,}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_R_GT_OBS'][nominal_alpha_idx]) + '\n'
        report += "\t\t[INFO] Original number of significant interactions: " \
                  + "{:,}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx]) + '\n'
        report += "\t\t[INFO] Mean number of significant randomized interactions: " \
                  + "{0:,.2f}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx]) + '\n'
        report += "\t\t[INFO] Standard deviation of significant randomized interactions: " \
                  + "{:.2f}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_R_SD'][nominal_alpha_idx]) + '\n'
        report += "\t\t[INFO] Z-score: "
        if self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx] != "NA":
            report += "{:.2f}".format(self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx]) + '\n'
        else:
            report += "NA (Increase number of iterations!)" + '\n'
        report += "\t\t[INFO] Estimated FDR: " \
                  + "{:.5f}".format(self._randomization_info_dict['RESULTS']['FDR'][nominal_alpha_idx]) + '\n'

        report += "[INFO] End of report." + '\n'

        return report

    def get_randomization_info_table_row(self, description: str = "DESCRIPTION", nominal_alphas_selected: [float] = None):
        """
        :return: String consisting of a header line and a line with values relating to last performed randomization
        analysis
        """

        # Truncate description
        description_truncated = (description[:30] + " ...") if len(description) > 30 else description

        # Get list of available nominal alphas
        nominal_alphas = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']

        # Check if there are any results for the passed list of nominal alphas
        if nominal_alphas_selected is not None:
            for nom_a in nominal_alphas_selected:
                if nom_a not in nominal_alphas:
                    report = "[ERROR] No results available for a nominal alpha of " + "{:.5f}".format(nom_a) + '!' + '\n'
                    report += "\t[ERROR] Results available for: " + ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas) +  '\n'
                    return report

        # Header row
        # ----------

        # Table tag and prefix for output
        table_row = ":TR_RANDOM:" + '\t'
        table_row += "DESCRIPTION" + '\t'

        # Input parameters
        table_row += "INPUT_I_NUM" + '\t'
        table_row += "ITER_NUM" + '\t'
        table_row += "RANDOM_SEED" + '\t'
        table_row += "NOMINAL_ALPHA" + '\t'

        # Results
        table_row += "SIG_NUM_R_GT_OBS" + '\t'
        table_row += "SIG_NUM_O" + '\t'
        table_row += "SIG_NUM_R_MEAN" + '\t'
        table_row += "SIG_NUM_R_SD" + '\t'
        table_row += "Z_SCORE" + '\t'
        table_row += "FDR" + '\n'

        # Row with values
        # ---------------

        for nominal_alpha in nominal_alphas:
            if nominal_alphas_selected is not None and nominal_alpha not in nominal_alphas_selected:
                continue
            nominal_alpha_idx = nominal_alphas.index(nominal_alpha)

            # Extract and format values from dictionary
            input_i_num = str(self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0])
            iter_num = str(self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0])
            random_seed = str(self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0])
            nominal_alpha = "{:.5f}".format(
                self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA'][nominal_alpha_idx])

            sig_num_r_gto = str(self._randomization_info_dict['RESULTS']['SIG_NUM_R_GT_OBS'][nominal_alpha_idx])
            sig_num_o = str(self._randomization_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_idx])
            sig_num_r_mean = "{:.2f}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_idx])
            sig_num_r_sd = "{:.2f}".format(self._randomization_info_dict['RESULTS']['SIG_NUM_R_SD'][nominal_alpha_idx])
            if self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx] != "NA":
                z_score = "{:.2f}".format(self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_idx])
            else:
                z_score = "NA"
            fdr = "{:.5f}".format(self._randomization_info_dict['RESULTS']['FDR'][nominal_alpha_idx])

            # Table tag and prefix for output
            table_row += ":TR_RANDOM:" + '\t'
            table_row += description_truncated + '\t'

            # Input parameters
            table_row += input_i_num + '\t'
            table_row += iter_num + '\t'
            table_row += random_seed + '\t'
            table_row += nominal_alpha + '\t'

            # Results
            table_row += sig_num_r_gto + '\t'
            table_row += sig_num_o + '\t'
            table_row += sig_num_r_mean + '\t'
            table_row += sig_num_r_sd + '\t'
            table_row += z_score + '\t'
            table_row += fdr + '\n'

        return table_row

    def _get_rp_inter_dict(self,
                           interaction_set: DiachromaticInteractionSet = None,
                           min_rp_num: int = 0):
        """
        Create dictionary required needed for randomization of interactions.

        :param min_rp_num: Interactions with fewer read pairs are not taken into account
        :return: Dictionary that contains the number of interactions (values) for all occurring read pair numbers
        (keys) for the interaction set of this object.
        """

        rp_inter_dict = {}
        for d_inter in interaction_set.interaction_list:
            if min_rp_num <= d_inter.rp_total:
                if d_inter.rp_total in rp_inter_dict:
                    rp_inter_dict[d_inter.rp_total] += 1
                else:
                    rp_inter_dict[d_inter.rp_total] = 1
        return rp_inter_dict

    def _get_list_of_p_values_from_randomized_data(self,
                                                   interaction_set: DiachromaticInteractionSet = None,
                                                   random_seed: int = None):
        """
        This function generates randomized simple and twisted read pair counts for all interactions of the
        interaction set of this object and calculates associated P-values.

        :return: List of P-values for randomized interactions
        """

        # Init random generator
        random.seed(random_seed)

        random_pval_list = []

        for rp_num, i_num in self._rp_inter_dict.items():

            # Generate random simple read pair counts for current read pair number
            simple_count_list = list(binom.rvs(rp_num, p=0.5, size=i_num))

            # Calculate P-values and append to list
            for simple_count in simple_count_list:
                nnl_pval = self.p_values.get_binomial_nnl_p_value(simple_count, rp_num - simple_count)
                random_pval_list.append(nnl_pval)

        return random_pval_list

    def _perform_one_iteration(self, nnl_nominal_alphas: [float] = None, random_seed: int = None):
        """
        This function performs a single iteration of the randomization procedure.

        :param nnl_nominal_alphas: List of nominal alphas (negative natural logarithm)
        :param random_seed: Number used to init random generator
        :return: List of numbers of randomized significant interactions for each nominal alpha
        """

        # Get list of P-values for randomized interactions
        randomized_nnl_pvals = self._get_list_of_p_values_from_randomized_data(random_seed=random_seed)

        # Determine number of randomized significant interactions for each nominal alpha in 'nnl_nominal_alphas'
        sig_num_r_list = self._determine_significant_pvals_at_nominal_alphas_nnl(nnl_nominal_alphas=nnl_nominal_alphas,
                                                                                 nnl_p_values=randomized_nnl_pvals)

        return sig_num_r_list

    def _perform_n_iterations(self, iter_start_idx: int, n, nnl_nominal_alphas: [float], verbose: bool=False):
        """
        This function performs a given number of iterations of the randomization procedure.

        :param iter_start_idx: First iteration index
        :param n: Number of iterations performed in the function call
        :param nnl_nominal_alphas: List of nominal alphas (negative natural logarithm)
        :param verbose: If true, messages about progress will be written to the screen
        :return: List of lists with numbers of randomized significant interactions. Each list contains the numbers of
        significant interactions for the different nominal alphas in 'nnl_nominal_alphas'.
        """

        # Iteration indices are added to the subordinate random seed
        iter_idx_range = range(iter_start_idx, n)

        # Init list of lists with numbers of randomized significant interactions for each iteration
        sig_num_r_list_of_lists = []

        if verbose:
            print("\t\t[INFO] Performing " + str(len(iter_idx_range)) + " iterations ...")
            print("\t\t\t[INFO] First iteration indices: " + ", ".join(str(i) for i in iter_idx_range[:10]) + ", ...")

        # Perform each iteration with its own random seed that corresponds by adding the iteration index
        for iter_idx in iter_idx_range:
            sig_num_r_list_of_lists.append(self._perform_one_iteration(nnl_nominal_alphas=nnl_nominal_alphas,
                                                                       random_seed=self._random_seed + iter_idx))

        return sig_num_r_list_of_lists

    def _determine_significant_pvals_at_nominal_alphas_nnl(self, nnl_nominal_alphas: [float], nnl_p_values: [float]):
        """
        Determine numbers of significant P-values at different nominal alphas (thresholds).

        :param nnl_nominal_alphas: List of nominal alphas (negative natural logarithm)
        :param nnl_p_values: List of P-values (negative natural logarithm)
        :return: List which has the same length as 'nominal_alphas' and contains, for each nominal alpha 'a',
        that number of P-values that are smaller or equal than 'a'.
        """

        # Sort input lists in ascending order
        nnl_nominal_alphas = sorted(nnl_nominal_alphas, reverse=True)
        nnl_p_values = sorted(nnl_p_values, reverse=True)

        # List of significant P-value numbers that will be returned
        sig_num_list = []

        # Index variable for the list of P-values
        nnl_p_value_idx = 0

        # Counting variable for the cumulative number of significant P-values
        sig_num = 0

        # Go through the list of nominal alphas sorted in ascending order
        for nnl_nominal_alpha in nnl_nominal_alphas:

            # As long as the P-values are smaller than the current nominal alpha,
            # increment index and counter variable
            while nnl_p_value_idx < len(nnl_p_values) and nnl_nominal_alpha <= nnl_p_values[nnl_p_value_idx]:
                sig_num += 1
                nnl_p_value_idx += 1

            # As  soon as the P-value is larger than the current nominal alpha,
            # add the current number of interactions to the list that will be returned
            sig_num_list.append(sig_num)

        return sig_num_list

    def get_randomization_info_plot(self,
                                    nominal_alpha_selected: float = None,
                                    pdf_file_name: str = "randomization.pdf",
                                    description: str = "DESCRIPTION"):
        """
        This function creates a graphical representation of the results from the last randomization analysis performed.

        :param nominal_alpha_selected: Nominal alpha for which the plot should be created
        :param pdf_file_name: Name of the PDF file that will be created
        :param description: Name of of the analysis that will be shown in the graphical representation
        :return: Nothing, if no FDR procedure has been performed yet or otherwise a 'Figure' object of 'matplotlib'
        that can be displayed in a Jupyter notebook
        """

        # Check whether an analysis has already been performed
        if self._randomization_info_dict is None:
            print("[ERROR] No analysis has been performed yet! There is nothing to plot.")
            return

        # Truncate description string
        description_truncated = (description[:30] + " ...") if len(description) > 30 else description

        # Check if there are results for the nominal alpha passed
        nominal_alphas_available = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']
        if nominal_alpha_selected not in nominal_alphas_available:
            if nominal_alpha_selected is None:
                nominal_alpha_selected_formatted = "None"
            else:
                nominal_alpha_selected_formatted = "{:.5f}".format(nominal_alpha_selected)
            report = "[ERROR] No results available for a nominal alpha of: " + nominal_alpha_selected_formatted + '!' + '\n'
            report += "\t[ERROR] Results available for: " + ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas_available) +  '\n'
            print(report)
            return

        # Set common plot parameters
        header_font_size = 10
        bin_face_color = 'blue'
        bin_edge_color = 'darkblue'

        # Extract data from randomization_info_dict
        # -----------------------------------------

        # Input parameters

        nominal_alphas = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']
        if len(nominal_alphas) < 5:
            nominal_alphas_formatted = ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas)
        else:
            nominal_alphas_formatted =  ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas[:3]) + ", ..., " + "{:.5f}".format(nominal_alphas[-1])
        iter_num = self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0]
        input_interactions_num = self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]
        random_seed = self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]

        # Results
        nominal_alpha_selected_idx = nominal_alphas.index(nominal_alpha_selected)
        if self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_selected_idx] == "NA":
            report = "[ERROR] There is no Z-Score available. Cannot create plot!" + '\n'
            print(report)
            return
        sig_num_r_list = self._randomization_info_dict['SIG_NUM_R_LISTS'][nominal_alpha_selected]
        sig_num_o = self._randomization_info_dict['RESULTS']['SIG_NUM_O'][nominal_alpha_selected_idx]
        i_num_randomized = self._randomization_info_dict['I_NUM_RANDOMIZED']
        sig_num_r_gt_obs = self._randomization_info_dict['RESULTS']['SIG_NUM_R_GT_OBS'][nominal_alpha_selected_idx]
        sig_num_r_mean = self._randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN'][nominal_alpha_selected_idx]
        sig_num_r_sd = self._randomization_info_dict['RESULTS']['SIG_NUM_R_SD'][nominal_alpha_selected_idx]
        z_score = float(self._randomization_info_dict['RESULTS']['Z_SCORE'][nominal_alpha_selected_idx])
        fdr = float(self._randomization_info_dict['RESULTS']['FDR'][nominal_alpha_selected_idx])

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
        ax[0].text(-0.2, 1.00, 'Randomization results', fontsize=18, fontweight='bold')
        ax[0].text(-0.18, 0.90, 'Analysis name: ' + description_truncated, fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.85, 'Number of input interactions: ' + "{:,}".format(input_interactions_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.80, 'Number of randomized interactions: ' + "{:,}".format(i_num_randomized),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.75, 'Number of iterations: ' + "{:,}".format(iter_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.70, 'Random seed: ' + str(random_seed),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.65, 'Available nominal alphas: ' + str(nominal_alphas_formatted),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.60, 'Selected nominal alpha: ' + "{:.5f}".format(nominal_alpha_selected),
                   fontsize=header_font_size, fontweight='bold')

        ax[0].text(-0.18, 0.50, 'Observed number of significant interactions: ' + "{:,}".format(sig_num_o),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.45,
                   'Iterations with more significant interactions than observed: ' + "{:,}".format(sig_num_r_gt_obs),
                   fontsize=header_font_size, color=sig_num_r_gt_obs_color, fontweight=sig_num_r_gt_obs_fontweight)
        ax[0].text(-0.18, 0.40,
                   'Mean number of significant randomized interactions: ' + "{:,.2f}".format(sig_num_r_mean),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.35,
                   'Standard deviation of significant randomized interactions: ' + "{:.2f}".format(sig_num_r_sd),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.30, 'Z-score: ' + "{:.2f}".format(z_score), fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.25, 'Estimated FDR : ' + "{:.5f}".format(fdr), fontsize=header_font_size, fontweight='bold')

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
        y_pos_lab = max(n) - 1 * (max(n) / 15)
        x_pos_lab = sig_num_r_mean + 6 * sig_num_r_sd
        ax[2].text(x_pos_lab, y_pos_lab, 'Observed: ' + "{:,}".format(sig_num_o),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})
        y_pos_lab = max(n) - 3 * (max(n) / 15)
        ax[2].text(x_pos_lab, y_pos_lab, 'Z-score: ' + "{:.2f}".format(z_score),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})
        y_pos_lab = max(n) - 5 * (max(n) / 15)
        ax[2].text(x_pos_lab, y_pos_lab, 'FDR: ' + "{:.5f}".format(fdr),
                   bbox={'color': 'w', 'alpha': 0.5, 'pad': 4})

        # Finalize save to PDF and return 'Figure' object
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    def get_randomization_info_plot_at_chosen_fdr_threshold(self,
                                                            pdf_file_name: str = "fdr_plot.pdf",
                                                            description: str = "DESCRIPTION",
                                                            chosen_fdr_threshold: float = 0.05):
        """
        This function creates a graphical representation of the results of the last randomization procedure performed.

        :param pdf_file_name: Name of the PDF file that will be created
        :param description: Name of of the analysis that will be shown in the graphical representation
        :param chosen_fdr_threshold: We are looking for the largest nominal alpha with an FDR below this threshold.
        :return: Nothing, if no randomization procedure has been performed yet or otherwise a 'Figure' object of
        'matplotlib' that can be displayed in a Jupyter notebook
        """

        # Check whether an analysis has already been performed
        if self._randomization_info_dict is None:
            print("[ERROR] No analysis has been performed yet! There is nothing to plot.")
            return

        # Set common plot parameters
        hv_lwd = 0.5  # line width of horizontal and vertical red dashed lines
        hv_col = 'red'  # color of horizontal and vertical red dashed lines
        header_font_size = 10

        # Truncate string for analysis name
        description_truncated = (description[:30] + " ...") if len(description) > 30 else description

        # Extract data from fdr_info_dict
        # -------------------------------

        # Input parameters
        chosen_fdr_thresh = chosen_fdr_threshold
        pval_thresh_max = max(self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA'])
        iter_num = self._randomization_info_dict['INPUT_PARAMETERS']['ITER_NUM'][0]
        random_seed = self._randomization_info_dict['INPUT_PARAMETERS']['RANDOM_SEED'][0]
        input_interactions_num = self._randomization_info_dict['INPUT_PARAMETERS']['INPUT_INTERACTIONS_NUM'][0]
        nominal_alphas = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']
        if len(nominal_alphas) < 5:
            nominal_alphas_formatted = ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas)
        else:
            nominal_alphas_formatted =  ", ".join("{:.5f}".format(nominal_alpha) for nominal_alpha in nominal_alphas[:3]) + ", ..., " + "{:.5f}".format(nominal_alphas[-1])

        i_num_randomized = self._randomization_info_dict['I_NUM_RANDOMIZED']

        # Determine results index for chosen FDR threshold
        result_index = 0
        for fdr_idx in range(0, len(self._randomization_info_dict['RESULTS']['FDR'])):
            if self._randomization_info_dict['RESULTS']['FDR'][fdr_idx] <= chosen_fdr_thresh:
                result_index = fdr_idx

        # Results
        pval_thresh_column = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA']
        pval_thresh_result = self._randomization_info_dict['RESULTS']['NOMINAL_ALPHA'][result_index]

        # Determine minimum number of read pairs required for significance and associated P-value for each nominal alpha
        min_rp_num_column =  []
        min_rp_num_pval_column = []
        for nom_alpha in pval_thresh_column:
            min_rp_num, min_rp_num_pval = self.p_values.find_smallest_significant_n(nom_alpha)
            min_rp_num_column.append(min_rp_num)
            min_rp_num_pval_column.append(min_rp_num_pval)
        min_rp_num_result = min_rp_num_column[result_index]
        min_rp_num_pval_result = min_rp_num_pval_column[result_index]

        sig_num_o_column = self._randomization_info_dict['RESULTS']['SIG_NUM_O']
        sig_num_o_result = self._randomization_info_dict['RESULTS']['SIG_NUM_O'][result_index]

        sig_num_r_gt_obs = self._randomization_info_dict['RESULTS']['SIG_NUM_R_GT_OBS'][result_index]

        # Highlight text if there are iterations with more significant interactions than originally observed
        if 0 < sig_num_r_gt_obs:
            sig_num_r_gt_obs_color = 'red'
            sig_num_r_gt_obs_fontweight = 'bold'
        else:
            sig_num_r_gt_obs_color = 'black'
            sig_num_r_gt_obs_fontweight = 'normal'

        sig_num_r_mean_column = self._randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN']
        sig_num_r_mean_result = self._randomization_info_dict['RESULTS']['SIG_NUM_R_MEAN'][result_index]

        sig_num_r_sd_result = self._randomization_info_dict['RESULTS']['SIG_NUM_R_SD'][result_index]

        z_score_column = self._randomization_info_dict['RESULTS']['Z_SCORE']
        z_score_result = self._randomization_info_dict['RESULTS']['Z_SCORE'][result_index]

        if "NA" in z_score_column:
            plot_z_scores = False
        else:
            plot_z_scores = True

        if z_score_result == "NA":
            z_score_result_formatted = "NA"
        else:
            z_score_result_formatted = "{0:,.2f}".format(z_score_result)

        fdr_column = self._randomization_info_dict['RESULTS']['FDR']
        fdr_result = self._randomization_info_dict['RESULTS']['FDR'][result_index]
        fdr_max = max(fdr_column)

        # Create figure with plots for individual columns
        # -----------------------------------------------

        fig, ax = plt.subplots(7, figsize=(7, 21.71), gridspec_kw={'height_ratios': [2, 1, 1, 1, 1, 1, 1]})

        # Add field with information about analysis
        plt.plot(ax=ax[0])
        ax[0].spines['left'].set_color('white')
        ax[0].spines['right'].set_color('white')
        ax[0].spines['top'].set_color('white')
        ax[0].spines['bottom'].set_color('white')
        ax[0].tick_params(axis='x', colors='white')
        ax[0].tick_params(axis='y', colors='white')
        ax[0].text(-0.2, 1.00, 'FDR results', fontsize=18, fontweight='bold')
        ax[0].text(-0.18, 0.90, 'Description: ' + description_truncated, fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.85, 'Number of input interactions: ' + "{:,}".format(input_interactions_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.80, 'Number of randomized interactions: ' + "{:,}".format(i_num_randomized),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.75, 'Number of iterations: ' + "{:,}".format(iter_num),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.70, 'Random seed: ' + str(random_seed),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.65, 'Available nominal alphas: ' + nominal_alphas_formatted, fontsize=header_font_size)
        ax[0].text(-0.18, 0.60, 'Chosen FDR threshold: ' + "{:.5f}".format(chosen_fdr_thresh),
                   fontsize=header_font_size, fontweight='bold')

        ax[0].text(-0.18, 0.50, 'Determined P-value threshold: ' + "{:.5f}".format(pval_thresh_result),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.14, 0.45, '-ln(' + "{:.5f}".format(pval_thresh_result) + ") = " + "{:.2f}".format(-log(pval_thresh_result)),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.40, 'Minimum read pair number: ' + str(min_rp_num_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.35,
                   'Smallest possible P-value with ' + str(min_rp_num_result) + ' read pairs: ' + "{:.5f}".format(
                       min_rp_num_pval_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.30, 'Observed number of significant interactions: ' + "{:,}".format(sig_num_o_result),
                   fontsize=header_font_size, fontweight='bold')
        ax[0].text(-0.18, 0.25,
                   'Iterations with more significant interactions than observed: ' + "{:,}".format(sig_num_r_gt_obs),
                   fontsize=header_font_size, color=sig_num_r_gt_obs_color, fontweight=sig_num_r_gt_obs_fontweight)
        ax[0].text(-0.18, 0.20, 'Mean number of randomized significant interactions: ' + "{0:,.2f}".format(sig_num_r_mean_result),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.15, 'Standard deviation of randomized significant interactions: ' + "{0:,.2f}".format(sig_num_r_sd_result),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.10, 'Z-score: ' + z_score_result_formatted,
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.05, 'Estimated FDR: ' + "{:.5f}".format(fdr_result), fontsize=header_font_size,
                   fontweight='bold')

        # Plot P-value thresholds
        ax[1].plot(pval_thresh_column, -log(pval_thresh_column))
        ax[1].set_title("P-value threshold: " + "{:.5f}".format(pval_thresh_result), loc='left')
        ax[1].set(xlabel="P-value threshold")
        ax[1].set(ylabel="-ln(P-value threshold)")
        ax[1].axhline(-log(pval_thresh_result), linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[1].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot minimum read pair numbers
        ax[2].plot(pval_thresh_column, min_rp_num_column)
        ax[2].set_title('Minimum read pair number: ' + str(min_rp_num_result), loc='left')
        ax[2].set(xlabel='P-value threshold')
        ax[2].set(ylabel='Read pair number')
        ax[2].axhline(min_rp_num_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[2].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Smallest possible P-values with minimum read pair numbers
        ax[3].plot(pval_thresh_column, min_rp_num_pval_column)
        ax[3].set_title('Smallest P-value with ' + str(min_rp_num_result) + ' read pairs: ' + "{:.5f}".format(
            min_rp_num_pval_result), loc='left')
        ax[3].set(xlabel='P-value threshold')
        ax[3].set(ylabel='Smallest P-value')
        ax[3].axhline(min_rp_num_pval_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[3].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot number of significant interactions
        ax[4].plot(pval_thresh_column, sig_num_r_mean_column, label='SIG_NUM_R')
        ax[4].plot(pval_thresh_column, sig_num_o_column, label='SIG_NUM_O')
        ax[4].set_title('Number of significant interactions: ' + "{:,}".format(sig_num_o_result) + ' (' + "{:,}".format(
            sig_num_r_mean_result) + ')', loc='left')
        ax[4].set(xlabel='P-value threshold')
        ax[4].set(ylabel='Significant interactions')
        ax[4].legend(loc="upper left", fontsize=8)
        ax[4].axhline(sig_num_o_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[4].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        # Plot Z-score
        if plot_z_scores:
            ax[5].plot(pval_thresh_column, z_score_column)
            ax[5].set_title('Z-score: ' + z_score_result_formatted, loc='left')
            ax[5].set(xlabel='P-value threshold')
            ax[5].set(ylabel='Z-score')
            ax[5].axhline(z_score_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
            ax[5].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        else:
            plt.plot(ax=ax[5])
            ax[5].spines['left'].set_color('white')
            ax[5].spines['right'].set_color('white')
            ax[5].spines['top'].set_color('white')
            ax[5].spines['bottom'].set_color('white')
            ax[5].tick_params(axis='x', colors='white')
            ax[5].tick_params(axis='y', colors='white')
            ax[5].text(-0.18, 0.90, 'The standard deviation of significant randomized interactions',
                       fontsize=header_font_size, fontweight='bold')
            ax[5].text(-0.18, 0.80, 'is zero for at least one nominal alpha!',
                       fontsize=header_font_size, fontweight='bold')
            ax[5].text(-0.18, 0.70, 'Therefore, the Z-scores cannot be plotted.',
                       fontsize=header_font_size)
            ax[5].text(-0.18, 0.60, 'The standard deviation is always zero, if only one iteration was performed.',
                       fontsize=header_font_size)
            ax[5].text(-0.18, 0.50, 'If only a few iterations have been performed,',
                       fontsize=header_font_size)
            ax[5].text(-0.18, 0.40, 'the standard deviation can sometimes also be zero.',
                       fontsize=header_font_size)
            ax[5].text(-0.18, 0.30, 'Try a larger number of iterations!.',
                       fontsize=header_font_size)

        # Plot estimated FDR
        ax[6].plot(pval_thresh_column, fdr_column)
        ax[6].set_title('Estimated FDR: ' + "{:.5f}".format(fdr_result), loc='left')
        ax[6].set(xlabel='P-value threshold')
        ax[6].set(ylabel='FDR')
        ax[6].axhline(fdr_result, linestyle='--', color=hv_col, linewidth=hv_lwd)
        ax[6].axvline(pval_thresh_result, linestyle='--', color=hv_col, linewidth=hv_lwd)

        ax[6].text(pval_thresh_max - pval_thresh_max / 5, fdr_result + fdr_max / 16,
                   'FDR: ' + "{:.5f}".format(fdr_result),
                   bbox={'color': 'lightblue', 'alpha': 0.5, 'pad': 4})

        ax[6].text(pval_thresh_result + pval_thresh_max / 60, fdr_result - fdr_max / 9,
                   'P-value: ' + "{:.5f}".format(pval_thresh_result),
                   bbox={'color': 'lightblue', 'alpha': 0.5, 'pad': 4})

        # Finalize save to PDF and return 'Figure' object
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig