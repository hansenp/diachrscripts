from scipy.stats import binom
from numpy import arange, log, random
import matplotlib.pyplot as plt
import warnings
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

warnings.filterwarnings('always')


class RandomizeInteractionSet:

    def __init__(self, random_seed: int = 42, interaction_set: DiachromaticInteractionSet = None):

        # Set random seed to be able to reproduce results
        self._random_seed = random_seed
        random.seed(random_seed)

        # Interaction set on which analyzes are performed
        self._interaction_set = interaction_set

        # Dictionary with all information on the last FDR procedure performed
        self._fdr_info_dict = None

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
        This function executes the FDR procedure.

        :param chosen_fdr_thresh: Predefined FDR threshold
        :param pval_thresh_max: Largest P-value threshold
        :param pval_thresh_step_size: Step size to increase the P-value threshold
        :param verbose: If true, messages about progress will be written to the screen
        :return:
        """

        if verbose:
            print("[INFO] Performing FDR procedure ...")

        # Set random seed
        random.seed(self._random_seed)

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
            nnnl_pval = self._interaction_set._p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)
            pvals_observed.append(nnnl_pval)

        if verbose:
            print("\t[INFO] Getting list of randomized P-values ...")

        # Get dictionary that stores the numbers of interactions with n read pairs and list of random P-values
        rp_inter_dict = self._get_rp_inter_dict()
        pvals_randomized = self._get_list_of_p_values_from_randomized_data(rp_inter_dict=rp_inter_dict)

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
                    'INPUT_INTERACTIONS_NUM': [len(self._interaction_set.interaction_list)]
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

    def perform_randomization_analysis(self, nominal_alpha: float = 0.01, n_iter: int = 1000):

        # Check input parameters
        if nominal_alpha <= 0 or 1.0 < nominal_alpha:
            raise ValueError("Nominal alpha must be in ]0,1]!")

        # Determine number of significant interactions at nominal alpha

        # Get dictionary that stores the numbers of interactions with n read pairs

        # Use dictionary to get list of P-values for randomized data and determine number of significant interactions

    def _get_rp_inter_dict(self):
        """
        :return: Dictionary that contains the number of interactions (values) for all occurring read pair numbers
        (keys) for the interaction set of this object.
        """

        rp_inter_dict = {}
        for d_inter in self._interaction_set.interaction_list:
            if d_inter.rp_total in rp_inter_dict:
                rp_inter_dict[d_inter.rp_total] += 1
            else:
                rp_inter_dict[d_inter.rp_total] = 1
        return rp_inter_dict

    def _get_list_of_p_values_from_randomized_data(self, rp_inter_dict: dict = None):
        """
        This function generates randomized simple and twisted read pair counts for all interactions of  the
        interaction set of this object and calculates associated P-values.

        :param rp_inter_dict: Dictionary generated with function '_get_rp_inter_dict'
        :return: List of P-values for randomized interactions
        """

        random_pval_list = []

        for rp_num, i_num in rp_inter_dict.items():

            # Generate random simple read pair counts for current read pair number
            simple_count_list = list(binom.rvs(rp_num, p=0.5, size=i_num))

            # Calculate P-values and append to list
            for simple_count in simple_count_list:
                pval = self._interaction_set._p_values.get_binomial_nnl_p_value(simple_count, rp_num - simple_count)
                random_pval_list.append(pval)

        return random_pval_list

    def get_fdr_info_plot(self, pdf_file_name: str = None, analysis_name: str = None):
        """
        This function creates a graphical representation of the results from the last FDR procedure performed.

        :param pdf_file_name: Name of the PDF file that will be created
        :param analysis_name: Name of of the analysis that will be shown in the graphical representation
        :return: Nothing, if no FDR procedure has been performed yet or otherwise a 'Figure' object of matplotlib
        that can be displayed in a Jupyter notebook
        """

        # Check whether an analysis has already been performed
        if self._fdr_info_dict is None:
            print("[ERROR] No analysis has been performed yet! There is nothing to plot.")
            return

        # Set parameters that all plots have in common
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

        ax[0].text(-0.18, 0.60, 'Determined P-value threshold: ' + "{:.5f}".format(pval_thresh_result),
                   fontsize=header_font_size,
                   fontweight='bold')
        ax[0].text(-0.18, 0.55, 'Determined -ln(P-value threshold): ' + "{:.5f}".format(nnl_pval_thresh_result),
                   fontsize=header_font_size)
        ax[0].text(-0.18, 0.50, 'Minimum read pair number: ' + str(min_rp_num_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.45,
                   'Smallest possible P-value with ' + str(min_rp_num_result) + ' read pairs: ' + "{:.5f}".format(
                       min_rp_num_pval_result), fontsize=header_font_size)
        ax[0].text(-0.18, 0.40,
                   'Number of interactions with ' + str(min_rp_num_result) + ' or more read pairs: ' + "{:,}".format(
                       min_rp_num_inter_num), fontsize=header_font_size, fontweight='bold', color=warn_min_rp_inter_col)
        ax[0].text(-0.18, 0.35, 'Number of significant interactions: ' + "{:,}".format(sig_num_o_result),
                   fontsize=header_font_size,
                   fontweight='bold')
        ax[0].text(-0.18, 0.30, 'Number of randomized significant interactions: ' + "{:,}".format(sig_num_r_result),
                   fontsize=header_font_size, color=warn_rand_sig_inter_col)
        ax[0].text(-0.18, 0.25, 'Estimated FDR: ' + "{:.5f}".format(fdr_result), fontsize=header_font_size,
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

        # Finalize save to PDF and return Figure object
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig
