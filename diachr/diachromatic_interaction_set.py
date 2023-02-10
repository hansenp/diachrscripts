import gzip
import os
import copy
import warnings
import random
from numpy import arange, log10
from collections import defaultdict
from .diachromatic_interaction import DiachromaticInteraction
from .diachromatic_interaction import DiachromaticInteraction11
from diachr.binomial_model import BinomialModel

warnings.filterwarnings('always')


class DiachromaticInteractionSet:
    """
    This class can read interactions from one ore more Diachromatic interaction files and perform the following
    operations on them:

        1. Calculation of P-values and classification into unbalanced and balanced
        2. Selection of two reference interaction sets, one for balanced and one for unbalanced
        3. Writing interactions that occur in a required number of input files to a new file

    If interactions have been evaluated and categorized (1,2), the output format is expanded by two columns on the left.
    Column 10 contains the negative of the natural logarithm of the P-value and column 11 the interaction category
    (DIX, DI, UI or UIR).
    """

    def __init__(self, enriched_digests_file: str = None, rpc_rule: str = 'ht'):

        # Dictionary that contains all interaction objects
        self._inter_dict = defaultdict(DiachromaticInteraction)

        # Rule of how to use the four read pair counts of interactions in this set
        if not (rpc_rule == 'st' or rpc_rule == 'ht'):
            raise ValueError(
                "Invalid tag for rule how to use read pair counts: " + rpc_rule + "! Use \'st\' or \'ht\'.")
        else:
            self.rpc_rule = rpc_rule

        # Object for calculating P-values
        if rpc_rule == 'st':
            self._p_values = BinomialModel(two_sided=True)
        else:
            self._p_values = BinomialModel(two_sided=False)

        # Smallest P-value threshold used so far
        self._smallest_pval_thresh = 1.0

        # Dictionary with information about read files
        self._read_file_info_dict = {'I_FILE': [], 'I_NUM': [], 'MIN_RP_NUM': [], 'MIN_DIST': [],
                                     'I_NUM_SKIPPED_RP': [], 'I_NUM_SKIPPED_DIST': [], 'I_NUM_ADDED': [],
                                     'I_SET_SIZE': []}

        # Dictionary with information about the last writing process
        self._write_file_info_dict = {}

        # Dictionary with information about the last evaluation and categorization process
        self._eval_cat_info_dict = {}

        # Dictionary with information about the last reference selection process
        self._select_ref_info_dict = {}

        # Take information about digests selected for enrichment from BED file
        # By default, the enrichment states of digests are taken from the interaction files.
        # If this file is passed, the enrichment sates are taken from this file.
        self._enriched_digests_set = None
        if enriched_digests_file is not None:
            print("[INFO] Reading list with digests selected for enrichment ...")
            self._enriched_digests_set = set()
            with open(enriched_digests_file, 'rt') as fp:
                for line in fp:
                    chr, sta, end = line.rstrip().split('\t')
                    self._enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))
            print("\t[INFO] Read " + "{:,}".format(len(self._enriched_digests_set)) + " digests ...")
            print("[INFO] ... done.")

    def parse_file(self, i_file: str = None, min_rp_num: int = 0, min_dist: int = 0, verbose: bool = False):
        """
        Parses a file with interactions. For interactions that have already been parsed from another file,
        only the read pair counts are added.

        :param i_file: Diachromatic interaction file
        :param min_rp_num: Interactions with a smaller number of read pairs are skipped
        :param min_dist: Interactions with shorter distances are skipped
        :param verbose: If true, messages about progress will be written to the screen
        """

        # Check if interaction file exists
        if not os.path.exists(i_file):
            raise ValueError("Could not find file at %s" % i_file)

        # Check if interaction file has already been read
        if i_file in self._read_file_info_dict['I_FILE']:
            warnings.warn("This interaction file has already been read!" + '\n' + \
                          "Filename: " + i_file + '\n' + \
                          "Won't add interactions from this file to the interaction set.")
            return

        if verbose:
            print("[INFO] Parsing Diachromatic interaction file ...")
            print("\t[INFO] " + i_file)

        if i_file.endswith(".gz"):
            with gzip.open(i_file, 'rt') as fp:
                n_lines = 0
                n_skipped_rp = 0
                n_skipped_dist = 0
                for line in fp:
                    n_lines += 1
                    if verbose:
                        if n_lines % 1000000 == 0:
                            print("\t[INFO] Parsed " + "{:,}".format(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.rp_total < min_rp_num:
                        n_skipped_rp += 1
                        continue
                    if d_inter.i_dist < min_dist:
                        n_skipped_dist += 1
                        continue
                    elif d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple_1=d_inter._simple_1,
                                                                              simple_2=d_inter._simple_2,
                                                                              twisted_1=d_inter._twisted_1,
                                                                              twisted_2=d_inter._twisted_2)
                    else:
                        self._inter_dict[d_inter.key] = d_inter
        else:
            with open(i_file) as fp:
                n_lines = 0
                n_skipped_rp = 0
                n_skipped_dist = 0
                for line in fp:
                    n_lines += 1
                    if verbose:
                        if n_lines % 1000000 == 0:
                            print("\t[INFO] Parsed " + "{:,}".format(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.rp_total < min_rp_num:
                        n_skipped_rp += 1
                        continue
                    if d_inter.i_dist < min_dist:
                        n_skipped_dist += 1
                        continue
                    elif d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple_1=d_inter._simple_1,
                                                                              simple_2=d_inter._simple_2,
                                                                              twisted_1=d_inter._twisted_1,
                                                                              twisted_2=d_inter._twisted_2)
                    else:
                        self._inter_dict[d_inter.key] = d_inter

        # Keep track of the name of the input file and the number of interactions
        self._read_file_info_dict['I_FILE'].append(i_file)
        self._read_file_info_dict['I_NUM'].append(n_lines)
        self._read_file_info_dict['MIN_RP_NUM'].append(min_rp_num)
        self._read_file_info_dict['MIN_DIST'].append(min_dist)
        self._read_file_info_dict['I_NUM_SKIPPED_RP'].append(n_skipped_rp)
        self._read_file_info_dict['I_NUM_SKIPPED_DIST'].append(n_skipped_dist)
        self._read_file_info_dict['I_NUM_ADDED'].append(n_lines - n_skipped_rp - n_skipped_dist)
        i_set_size = len(self._inter_dict)
        self._read_file_info_dict['I_SET_SIZE'].append(i_set_size)

        if verbose:
            print("\t[INFO] Set size: " + "{:,}".format(i_set_size))
            print("[INFO] ... done.")

    def _parse_line(self, line: str = None) -> DiachromaticInteraction:
        """
        Parses a Diachromatic interaction formatted line with 9 or 11 fields.

        :param line: Tab separated string 9 or 11 fields.
        :return: DiachromaticInteraction object
        """

        # Split string and check number of fields
        F = line.rstrip().split('\t')
        if len(F) < 9:
            raise ValueError("Malformed line {}".format(line))

        # Assign fields to variables
        chrA = F[0]
        fromA = int(F[1])
        toA = int(F[2])
        statusA = F[3]
        chrB = F[4]
        fromB = int(F[5])
        toB = int(F[6])
        statusB = F[7]
        st_string = F[8]  # something like 2:0:1:0, representing S1:S2:T1:T2, simple and twisted counts
        st_fields = st_string.split(":")
        if len(st_fields) != 4:
            raise ValueError("Malformed simple_1:simple_2:twisted_1:twisted_2 field in diachromatic line: " + line)
        else:
            rp_counts = [int(i) for i in st_fields]
            simple_1 = rp_counts[0]
            simple_2 = rp_counts[1]
            twisted_1 = rp_counts[2]
            twisted_2 = rp_counts[3]

        # Get enrichment status of digests from file for enriched digests if available
        if self._enriched_digests_set is not None:
            coord_key_da = chrA + "\t" + str(fromA) + "\t" + str(toA)
            coord_key_db = chrB + "\t" + str(fromB) + "\t" + str(toB)
            if coord_key_da in self._enriched_digests_set:
                statusA = 'E'
            else:
                statusA = 'N'

            if coord_key_db in self._enriched_digests_set:
                statusB = 'E'
            else:
                statusB = 'N'

        # The interaction read has not yet been evaluated and categorized
        if len(F) < 11:
            di_inter = DiachromaticInteraction(
                chrA=chrA,
                fromA=fromA,
                toA=toA,
                statusA=statusA,
                chrB=chrB,
                fromB=fromB,
                toB=toB,
                statusB=statusB,
                simple_1=simple_1,
                simple_2=simple_2,
                twisted_1=twisted_1,
                twisted_2=twisted_2)
            return di_inter

        # The interaction read has already been evaluated and categorized
        if len(F) == 11:
            log10_p_val = float(F[9])
            i_cat = F[10]
            di_11_inter = DiachromaticInteraction11(
                chrA=chrA,
                fromA=fromA,
                toA=toA,
                statusA=statusA,
                chrB=chrB,
                fromB=fromB,
                toB=toB,
                statusB=statusB,
                simple_1=simple_1,
                simple_2=simple_2,
                twisted_1=twisted_1,
                twisted_2=twisted_2,
                log10_pval=log10_p_val)

            # Set interaction category
            di_11_inter.set_category(i_cat)
            return di_11_inter

    def evaluate_and_categorize_interactions(self, pval_thresh: float = None, verbose: bool = False):
        """
        Calculate the P-value and define interaction category ('DI' or 'UI') for all interactions in this object.
        'DiachromaticInteraction' objects will be replaced by 'DiachromaticInteraction11' objects.

        :param pval_thresh: Interactions with a lower P-value will be classified as significant
        :param verbose:  If true, messages about progress will be written to the screen
        :return: Dictionary with information on processing.
        """

        # Check whether an equal or smaller P-value threshold was used before
        if self._smallest_pval_thresh <= pval_thresh:
            warnings.warn(
                "This interaction set has previously categorized using an equal or smaller threshold value!" + '\n'
                                                                                                               "Nothing is done. Interaction set remains unchanged.")
            return self._eval_cat_info_dict
        else:
            self._smallest_pval_thresh = pval_thresh

        if verbose:
            print("[INFO] Evaluate and categorize interactions ...")

        # Get smallest number of read pairs required for significance
        min_rp, min_rp_pval = self._p_values.find_smallest_significant_n(pval_thresh, verbose=False)

        # Transform P-value threshold
        log10_pval_thresh = -log10(pval_thresh)

        # Iterate interaction objects
        d11_inter_dict = {}
        n_di = 0
        n_ui = 0
        n_discarded = 0
        n_processed = 0
        for d_inter in self._inter_dict.values():
            n_processed += 1
            if verbose:
                if n_processed % 1000000 == 0:
                    print("\t[INFO] Processed " + "{:,}".format(n_processed) + " interactions ...")

            # Discard interactions that don't have enough read pairs to be significant
            if d_inter.rp_total < min_rp:
                n_discarded += 1
                continue

            # Get P-value
            if self.rpc_rule == 'st':
                log10_p_val = self._p_values.get_binomial_log10_p_value(d_inter.n_simple, d_inter.n_twisted)
            else:
                log10_p_val = self._p_values.get_binomial_log10_p_value(d_inter.n_heaviest_two, d_inter.n_lightest_two)

            # Determine interaction category
            if log10_pval_thresh <= log10_p_val:
                i_category = "DI"
                n_di += 1
            else:
                i_category = "UI"
                n_ui += 1

            # Create new 'DiachromaticInteraction11' with P-value and category
            di_11_inter = DiachromaticInteraction11(
                chrA=d_inter.chrA,
                fromA=d_inter._fromA,
                toA=d_inter._toA,
                statusA=d_inter.enrichment_status_tag_pair[0],
                chrB=d_inter.chrB,
                fromB=d_inter._fromB,
                toB=d_inter._toB,
                statusB=d_inter.enrichment_status_tag_pair[1],
                simple_1=d_inter._simple_1,
                simple_2=d_inter._simple_2,
                twisted_1=d_inter._twisted_1,
                twisted_2=d_inter._twisted_2,
                log10_pval=log10_p_val)

            # Set interaction category
            di_11_inter.set_category(i_category)
            d11_inter_dict[d_inter.key] = di_11_inter

        # Replace old list with new list of interactions
        self._inter_dict = d11_inter_dict

        # Prepare dictionary for report
        report_dict = {
            'PVAL_THRESH': [pval_thresh],
            'MIN_RP': [min_rp],
            'MIN_RP_PVAL': [min_rp_pval],
            'N_PROCESSED': [n_processed],
            'N_DISCARDED': [n_discarded],
            'N_BALANCED': [n_ui],
            'N_UNBALANCED': [n_di]
        }

        if verbose:
            print("[INFO] ... done.")

        self._eval_cat_info_dict = report_dict
        return report_dict

    def select_reference_interactions(self, selection_rule: str = "RPNUM", verbose: bool = False):
        """
        Select balanced reference interactions that match unbalanced interactions in terms of enrichment category and total number
        of read pairs per interaction and return a dictionary with information on this selection process.

        The difference to the original version is that it is not distinguished between NE and EN.
        Furthermore, an additional category DIX for DI interactions without counterpart in reference interactions
        is introduced.

        :return: Dictionary with interaction numbers in different categories
        """

        if verbose:
            print("[INFO] Select reference interactions ...")
            print("\t[INFO] Treating NE and EN as one category ...")

        if selection_rule not in ['RPNUM', 'RPMAX', 'RPMAX2', 'RPMAX3']:
            print('[ERROR] The selection rule must be \'RPNUM\', \'RPMAX\', \'RPMAX2\', \'RPMAX3\'!')
            return

        # Nested dictionary that stores the numbers of interactions (value) for different read pair numbers (key)
        rp_inter_dict = {'NN': {},
                         'NEEN': {},
                         'EE': {}}

        # Count number of interactions in different categories
        report_dict = {'NN': {'DIX': [0], 'DI': [0], 'UIR': [0], 'UI': [0]},
                       'NE': {'DIX': [0], 'DI': [0], 'UIR': [0], 'UI': [0]},
                       'EN': {'DIX': [0], 'DI': [0], 'UIR': [0], 'UI': [0]},
                       'EE': {'DIX': [0], 'DI': [0], 'UIR': [0], 'UI': [0]}}

        if verbose:
            print("\t[INFO] First pass: Count unbalanced interactions for different read pair counts ...")
        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() == 'DI':

                # Get enrichment status tag pair and read pair number
                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                report_dict[enrichment_pair_tag]['DI'][0] += 1
                if enrichment_pair_tag == 'NE' or enrichment_pair_tag == 'EN':
                    enrichment_pair_tag = 'NEEN'

                if selection_rule == 'RPNUM':
                    rp_total = d11_inter.rp_total
                elif selection_rule == 'RPMAX':
                    rp_total = max([d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX2':
                    rp_total = max([d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX3':
                    rpcs = sorted(
                        [d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2],
                        reverse=True)
                    rp_total = sum(rpcs[:2])
                else:
                    print('[ERROR] The selection rule must be \'RPNUM\', \'RPMAX\', \'RPMAX2\', \'RPMAX3\'!')
                    return

                if rp_total not in rp_inter_dict[enrichment_pair_tag]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] = 1
                else:
                    rp_inter_dict[enrichment_pair_tag][rp_total] += 1

        rp_inter_dict_before = copy.deepcopy(rp_inter_dict)
        ui_inter_dict = {'NN': 0,
                         'NEEN': 0,
                         'EE': 0}

        if verbose:
            print("\t[INFO] Second pass: Select balanced reference interactions for different read pair counts ...")

        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() != 'DI':

                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                if enrichment_pair_tag == 'NE' or enrichment_pair_tag == 'EN':
                    enrichment_pair_tag = 'NEEN'

                if selection_rule == 'RPNUM':
                    rp_total = d11_inter.rp_total
                elif selection_rule == 'RPMAX':
                    rp_total = max([d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX2':
                    rp_total = max([d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX3':
                    rpcs = sorted(
                        [d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2],
                        reverse=True)
                    rp_total = sum(rpcs[:2])
                else:
                    print('[ERROR] The selection rule must be \'RPNUM\', \'RPMAX\', \'RPMAX2\', \'RPMAX3\'!')
                    return

                if rp_total in rp_inter_dict[enrichment_pair_tag] and 0 < rp_inter_dict[enrichment_pair_tag][rp_total]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] -= 1
                    d11_inter.set_category('UIR')
                    report_dict[d11_inter.enrichment_status_tag_pair]['UIR'][0] += 1
                else:
                    ui_inter_dict[enrichment_pair_tag] += 1
                    d11_inter.set_category('UI')
                    report_dict[d11_inter.enrichment_status_tag_pair]['UI'][0] += 1

        if verbose:
            print("\t[INFO] Third pass: Moving DI interactions for which there is no reference to DIX ...")

        for d11_inter in self._inter_dict.values():
            if d11_inter.get_category() == 'DI':

                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                if enrichment_pair_tag == 'NE' or enrichment_pair_tag == 'EN':
                    enrichment_pair_tag = 'NEEN'

                if selection_rule == 'RPNUM':
                    rp_total = d11_inter.rp_total
                elif selection_rule == 'RPMAX':
                    rp_total = max([d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX2':
                    rp_total = max([d11_inter._twisted_1, d11_inter._twisted_2])
                elif selection_rule == 'RPMAX3':
                    rpcs = sorted(
                        [d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2],
                        reverse=True)
                    rp_total = sum(rpcs[:2])
                else:
                    print('[ERROR] The selection rule must be \'RPNUM\', \'RPMAX\', \'RPMAX2\', \'RPMAX3\'!')
                    return

                if rp_total in rp_inter_dict[enrichment_pair_tag] and 0 < rp_inter_dict[enrichment_pair_tag][rp_total]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] -= 1
                    d11_inter.set_category('DIX')
                    report_dict[d11_inter.enrichment_status_tag_pair]['DI'][0] -= 1
                    report_dict[d11_inter.enrichment_status_tag_pair]['DIX'][0] += 1

        # Sanity check DI and UIR must have the same number of interactions in each enrichment category
        for enr_cat in ['NN', 'EE']:
            if report_dict[enr_cat]['DI'][0] - report_dict[enr_cat]['UIR'][0] != 0:
                print("[ERROR] DI and UIR must have the same number of interactions in category NN and EE!")
                print("\t[ERROR] " + enr_cat)
                print("\t[ERROR] DI: " + str(report_dict[enr_cat]['DI'][0]))
                print("\t[ERROR] UIR: " + str(report_dict[enr_cat]['UIR'][0]))
                exit(1)
        neen_sum_di = report_dict['NE']['DI'][0] + report_dict['EN']['DI'][0]
        neen_sum_uir = report_dict['NE']['UIR'][0] + report_dict['EN']['UIR'][0]
        if neen_sum_di - neen_sum_uir != 0:
            print("[ERROR] DI and UIR must have the same number of interactions in category NE and EN!")
            print("\t[ERROR] DI-NE: " + str(report_dict['NE']['DI'][0] + report_dict['NE']['DI'][0]))
            print("\t[ERROR] DI-EN: " + str(report_dict['NE']['DI'][0] + report_dict['EN']['DI'][0]))
            print("\t[ERROR] DI-NEEN: " + str(neen_sum_di))
            print("\t[ERROR] UIR-NE: " + str(report_dict['NE']['UIR'][0] + report_dict['NE']['UIR'][0]))
            print("\t[ERROR] UIR-EN: " + str(report_dict['NE']['UIR'][0] + report_dict['EN']['UIR'][0]))
            print("\t[ERROR] UIR-NEEN: " + str(neen_sum_uir))
            return

        if verbose:
            print("[INFO] ... done.")

        self._select_ref_info_dict = report_dict
        return report_dict

    def deselect_reference_interactions(self):
        """
        This function undoes the selection of reference interactions. The tag 'DIX' is replaced with 'DI' and the tag
        'UIR' with 'UI'.
        """
        for d11_inter in self._inter_dict.values():
            if d11_inter.get_category() == 'DIX':
                d11_inter.set_category('DI')
            if d11_inter.get_category() == 'UIR':
                d11_inter.set_category('UI')

    def write_diachromatic_interaction_file(self, required_replicates: int = 1, target_file: str = None,
                                            verbose: bool = False):
        """
        Write interactions that occur in a specified minimum number of replicates to a file and return a dictionary
        with information on this writing process.

        :param required_replicates: Minimum number of replicates
        :param target_file_name: Generated file with interactions that occur in the required number of replicates
        :return: Dictionary with information on this writing process
        """

        if verbose:
            print("[INFO] Writing Diachromatic interaction file ...")
            print("\t[INFO] Required replicates: " + str(required_replicates))
            print("\t[INFO] Target file: " + target_file)

        # Write file
        out_fh = gzip.open(target_file, 'wt')
        n_complete_data = 0
        n_incomplete_data = 0
        for d_inter in self._inter_dict.values():
            if d_inter.has_data_for_required_replicate_num(required_replicates):
                n_complete_data += 1
                out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
            else:
                n_incomplete_data += 1
        out_fh.close()

        # Prepare dictionary for report
        report_dict = {
            'TARGET_FILE': [target_file],
            'REQUIRED_REPLICATES': [required_replicates],
            'N_INCOMPLETE_DATA': [n_incomplete_data],
            'N_COMPLETE_DATA': [n_complete_data]
        }

        if verbose:
            print("[INFO] ... done.")

        self._write_file_info_dict = report_dict
        return report_dict

    def shuffle_inter_dict(self, random_seed: int = None, verbose: bool = False):
        """
        Input files can be affected by sorting artifacts. This function shuffles the key value pairs and thus eliminates
        such artifacts.

        :param random_seed: Number used to init random generator
        :param verbose: If true, messages about progress will be written to the screen
        """

        if verbose:
            print("[INFO] Shuffling dictionary with interactions ...")

        if random_seed is not None:
            random.seed(random_seed)

        if verbose:
            print("\t[INFO] Random seed: " + str(random_seed))

        items = list(self._inter_dict.items())
        random.shuffle(items)
        self._inter_dict = defaultdict(DiachromaticInteraction, items)

        if verbose:
            print("[INFO] ... done.")

    def write_diachromatic_interaction_fdr_test_file(self, target_file: str = None, p_value_max: float = 0.05,
                                                     p_value_step: float = 0.00025, i_per_interval_requested: int = 10,
                                                     verbose: bool = False):
        """
        Generates a test file that contains the same number of interactions for consecutive P-value intervals.
        that can be used to test to check whether the P-value thresholds are used correctly during the FDR analysis.

        :param target_file: Output interaction file
        :param p_value_max: Maximum P-value
        :param p_value_step: P-value step size
        :param i_per_interval_requested: Requested number of interactions per interval
        :param verbose: If true, messages about progress will be written to the screen
        :return: Dictionary that contains intermediate results and summary statistics
        """

        if verbose:
            print("[INFO] Generating test file for FDR procedure ...")

        # Open interaction file for writing
        out_fh = open(target_file, 'wt')

        # Create list of ascending P-value thresholds
        p_threshs = arange(p_value_step, p_value_max + p_value_step, p_value_step)

        # Prepare lists for intermediate results
        pval_interval_list = []
        i_selected_list = []
        i_selected_cum_list = []

        # Iterate interaction set and write requested number of interactions for each P-value interval to file
        if verbose:
            print("\t[INFO] Iterating interaction set and writing interactions to test file ...")
        i_selected_cum = 0
        for p_thresh in p_threshs:
            i_selected = 0
            for d_inter in self.interaction_list:
                pval = self._p_values.get_binomial_p_value(d_inter.n_simple, d_inter.n_twisted)
                if (p_thresh - p_value_step) < pval and pval <= p_thresh:
                    i_selected += 1
                    out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
                    i_selected_cum += 1
                if i_selected == i_per_interval_requested:
                    break

            # Collect intermediate results
            pval_interval = ']' + "{:.5f}".format(p_thresh - p_value_step) + ';' + "{:.5f}".format(p_thresh) + ']'
            pval_interval_list.append(pval_interval)
            i_selected_list.append(i_selected)
            i_selected_cum_list.append(i_selected_cum)

            # If the requested number of interactions could not be selected for the current interval, issue a warning
            if i_selected < i_per_interval_requested:
                warnings.warn("[WARNING] Could not select the required number of interactions (only "
                              + str(i_selected) + " of " + str(i_per_interval_requested) +
                              ") for the P-value range " + pval_interval + '!')

        # Close interaction file
        out_fh.close()

        if verbose:
            print("[INFO] ... done.")

        # Prepare dictionary for report
        report_dict = {
            'TARGET_FILE': target_file,
            'RESULTS_TABLE':
                {'PVAL_INTERVAL': pval_interval_list,
                 'I_SELECTED': i_selected_list,
                 'I_SELECTED_CUM': i_selected_cum_list
                 }
        }
        return report_dict

    def get_read_file_info_dict(self):
        """
        :return: Dictionary that contains information about parsed interaction data.
        """

        return self._read_file_info_dict

    def get_read_file_info_report(self):
        """
        :return: String that contains information about parsed interaction data.
        """

        file_num = len(self._read_file_info_dict['I_FILE'])
        report = "[INFO] Report on reading files:" + '\n'
        report += "\t[INFO] Read interaction data from " + str(file_num) + " files:" + '\n'
        for i in range(0, file_num):
            report += "\t\t[INFO] " + "{:,}".format(self._read_file_info_dict['I_NUM'][i]) + " interactions from: \n" + \
                      "\t\t\t[INFO] " + self._read_file_info_dict['I_FILE'][i] + '\n' \
                                                                                 "\t\t\t[INFO] Minimum number of read pairs: " + str(
                self._read_file_info_dict['MIN_RP_NUM'][i]) + '\n' \
                                                              '\t\t\t[INFO] Skipped because less than ' + str(
                self._read_file_info_dict['MIN_RP_NUM'][i]) + \
                      ' read pairs: ' + \
                      "{:,}".format(self._read_file_info_dict['I_NUM_SKIPPED_RP'][i]) + '\n' + \
                      "\t\t\t[INFO] Minimum interaction distance: " + "{:,}".format(
                self._read_file_info_dict['MIN_DIST'][i]) + '\n' \
                                                            '\t\t\t[INFO] Skipped because shorter than ' + "{:,}".format(
                self._read_file_info_dict['MIN_DIST'][i]) + \
                      ' bp: ' + \
                      "{:,}".format(self._read_file_info_dict['I_NUM_SKIPPED_DIST'][i]) + '\n' + \
                      '\t\t\t[INFO] Added to set: ' + "{:,}".format(self._read_file_info_dict['I_NUM_ADDED'][i]) + '\n' \
                                                                                                                   '\t\t\t[INFO] Set size: ' + "{:,}".format(
                self._read_file_info_dict['I_SET_SIZE'][i]) + '\n'
        report += "\t[INFO] The interaction set has " + \
                  "{:,}".format(self._read_file_info_dict['I_SET_SIZE'][-1]) + " interactions." + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_read_file_info_table_row(self):
        """
        :return: String consisting of a header line and a lines with values relating to parsed interactions
        """

        # Header row
        table_row = ":TR_READ:" + "\t" + \
                    "I_FILE" + "\t" + \
                    "I_NUM" + "\t" + \
                    "MIN_RP_NUM" + "\t" + \
                    "I_NUM_SKIPPED_RP" + "\t" + \
                    "MIN_DIST" + "\t" + \
                    "I_NUM_SKIPPED_DIST" + "\t" + \
                    "I_NUM_ADDED" + "\t" + \
                    "I_SET_SIZE" + '\n'

        # Rows with values
        file_num = len(self._read_file_info_dict['I_FILE'])
        for i in range(0, file_num):
            table_row += ":TR_READ:" + "\t"
            table_row += self._read_file_info_dict['I_FILE'][i] + '\t'
            table_row += str(self._read_file_info_dict['I_NUM'][i]) + '\t'
            table_row += str(self._read_file_info_dict['MIN_RP_NUM'][i]) + '\t'
            table_row += str(self._read_file_info_dict['I_NUM_SKIPPED_RP'][i]) + '\t'
            table_row += str(self._read_file_info_dict['MIN_DIST'][i]) + '\t'
            table_row += str(self._read_file_info_dict['I_NUM_SKIPPED_DIST'][i]) + '\t'
            table_row += str(self._read_file_info_dict['I_NUM_ADDED'][i]) + '\t'
            table_row += str(self._read_file_info_dict['I_SET_SIZE'][i]) + '\n'

        return table_row

    def get_write_file_info_report(self):
        """
        :return: String that contains information about the last writing process.
        """

        report = "[INFO] Report on writing files:" + '\n'
        report += "\t[INFO] Wrote interactions that occur in at least " + \
                  str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + " replicates to:" + '\n'
        report += "\t\t[INFO] " + self._write_file_info_dict['TARGET_FILE'][0] + '\n'
        report += "\t[INFO] Interactions that occur in at least " + \
                  str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + " replicates: " + \
                  "{:,}".format(self._write_file_info_dict['N_COMPLETE_DATA'][0]) + '\n'
        report += "\t[INFO] Other interactions: " + \
                  "{:,}".format(self._write_file_info_dict['N_INCOMPLETE_DATA'][0]) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_write_file_info_table_row(self):
        """
        :return: String consisting of a header line and a line with values relating to written interactions
        """

        table_row = ":TR_WRITE:" + "\t" \
                                   "TARGET_FILE" + "\t" + \
                    "REQUIRED_REPLICATES" + "\t" + \
                    "N_INCOMPLETE_DATA" + "\t" + \
                    "N_COMPLETE_DATA" + '\n'
        table_row += ":TR_WRITE:" + "\t" + \
                     str(self._write_file_info_dict['TARGET_FILE'][0]) + "\t" + \
                     str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + "\t" + \
                     str(self._write_file_info_dict['N_INCOMPLETE_DATA'][0]) + "\t" + \
                     str(self._write_file_info_dict['N_COMPLETE_DATA'][0]) + '\n'

        return table_row

    def get_eval_cat_info_report(self):
        """
        :return: String that contains information about the last evaluation and categorization process.
        """

        report = "[INFO] Report on evaluation and categorization interactions:" + '\n'
        report += "\t[INFO] P-value threshold: " + "{:.5f}".format(
            self._eval_cat_info_dict['PVAL_THRESH'][0]) + '\n'
        report += "\t[INFO] Minimum number of read pairs required for significance: " + str(
            self._eval_cat_info_dict['MIN_RP'][0]) + '\n'
        report += "\t[INFO] Smallest P-value with " + str(self._eval_cat_info_dict['MIN_RP'][0]) + " read pairs: " + \
                  "{:.5f}".format(self._eval_cat_info_dict['MIN_RP_PVAL'][0]) + '\n'
        report += "\t[INFO] Processed interactions: " + "{:,}".format(self._eval_cat_info_dict['N_PROCESSED'][0]) + '\n'
        report += "\t[INFO] Discarded interactions: " + "{:,}".format(self._eval_cat_info_dict['N_DISCARDED'][0]) + '\n'
        report += "\t[INFO] Not significant interactions (UI): " + "{:,}".format(
            self._eval_cat_info_dict['N_BALANCED'][0]) + '\n'
        report += "\t[INFO] Significant interactions (DI): " + "{:,}".format(
            self._eval_cat_info_dict['N_UNBALANCED'][0]) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_eval_cat_info_table_row(self, out_prefix: str = None):
        """
        :return: String consisting of a header line and a line with values relating to evaluation and categorization
        of interactions
        """

        table_row = ":TR_EVAL_CAT:" + '\t' + \
                    "DESCRIPTION" + '\t' + \
                    "PVAL_THRESH" + '\t' + \
                    "MIN_RP" + '\t' + \
                    "MIN_RP_PVAL" + '\t' + \
                    "N_PROCESSED" + '\t' + \
                    "N_DISCARDED" + '\t' + \
                    "N_BALANCED" + '\t' \
                    "N_UNBALANCED" + '\n'

        table_row += ":TR_EVAL_CAT:" + '\t' + \
                     str(out_prefix) + '\t' + \
                     "{:.5f}".format(self._eval_cat_info_dict['PVAL_THRESH'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['MIN_RP'][0]) + '\t' + \
                     "{:.5f}".format(self._eval_cat_info_dict['MIN_RP_PVAL'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_PROCESSED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_DISCARDED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_BALANCED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_UNBALANCED'][0]) + '\n'

        return table_row

    def get_select_ref_info_dict(self):
        """
        :return: Dictionary that contains information about the last reference selection process.
        """
        return self._select_ref_info_dict

    def get_select_ref_info_report(self):
        """
        :return: String that contains information about the last reference selection process.
        """

        report = "[INFO] Report on selection of balanced reference interactions:" + '\n'
        report += "\t[INFO] Numbers of unbalanced interactions without balanced reference (DIX)" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict[enr_cat]['DIX'][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + "{:,}".format(
                self._select_ref_info_dict[enr_cat]['DIX'][0]) + '\n'
        report += "\t\t[INFO] Total: " + "{:,}".format(total) + '\n'
        report += "\t[INFO] Numbers of unbalanced interactions with balanced reference (DI)" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict[enr_cat]['DI'][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + "{:,}".format(
                self._select_ref_info_dict[enr_cat]['DI'][0]) + '\n'
        report += "\t\t[INFO] Total: " + "{:,}".format(total) + '\n'
        report += "\t[INFO] Numbers of balanced reference interactions (UIR)" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict[enr_cat]['UIR'][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + "{:,}".format(
                self._select_ref_info_dict[enr_cat]['UIR'][0]) + '\n'
        report += "\t\t[INFO] Total: " + "{:,}".format(total) + '\n'
        report += "\t[INFO] Numbers of balanced interactions not selected as reference (UI)" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict[enr_cat]['UI'][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + "{:,}".format(
                self._select_ref_info_dict[enr_cat]['UI'][0]) + '\n'
        report += "\t\t[INFO] Total: " + "{:,}".format(total) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_select_ref_info_table_row(self, description: str = None):
        """
        :return: String consisting of a header line and a line with values relating to the last selection of reference
        interactions
        """

        # Header line
        table_row = ":TR_SELECT:" + '\t'
        table_row += "DESCRIPTION" + '\t'
        for i_cat in ['DIX', 'DI', 'UIR', 'UI']:
            for enr_cat in ['NN', 'NE', 'EN', 'EE']:
                table_row += i_cat + '_' + enr_cat + '\t'
            if i_cat != 'UI':
                table_row += i_cat + '_TOTAL' + '\t'
            else:
                table_row += i_cat + '_TOTAL' + '\n'

        # Line with values
        table_row += ":TR_SELECT:" + '\t'
        table_row += str(description) + '\t'
        for i_cat in ['DIX', 'DI', 'UIR', 'UI']:
            total = 0
            for enr_cat in ['NN', 'NE', 'EN', 'EE']:
                total += self._select_ref_info_dict[enr_cat][i_cat][0]
                table_row += str(self._select_ref_info_dict[enr_cat][i_cat][0]) + '\t'
            if i_cat != 'UI':
                table_row += str(total) + '\t'
            else:
                table_row += str(total) + '\n'

        return table_row

    @property
    def interaction_list(self):
        return self._inter_dict.values()

    @property
    def interaction_dict(self):
        return self._inter_dict

    def get_num_of_inter_with_as_many_or_more_read_pairs(self, min_rp_num: int = None):
        """
        Determine the number of interactions that have a given number of more read pairs.

        :param min_rp_num: Minimum number of read pairs
        :return: Number of interactions with given number or more read pairs
        """
        i_num = 0
        for d_inter in self.interaction_list:
            if min_rp_num <= d_inter.rp_total:
                i_num += 1
        return i_num

    def sort_and_select_top_n_interactions(self, sort_by: str = 'RPNUM', top_n: int = 100000):

        # Get list of interactions
        interaction_list = list(self._inter_dict.values())

        if len(interaction_list) < top_n:
            print('[ERROR] Invalid \'top_n\' parameter! Must be less then the total number of interactions, ' + \
                  'which is ' + '{:,}'.format(len(interaction_list)) + '.')

        # Get list of interaction indices
        idx_list = list(range(0, len(interaction_list)))

        # Get RPNUM or RPMAX list
        if sort_by == 'RPNUM':
            rp_list = [d_inter.rp_total for d_inter in interaction_list]
        elif sort_by == 'RPMAX':
            rp_list = [max(d_inter._simple_1, d_inter._simple_2, d_inter._twisted_1, d_inter._twisted_2) \
                       for d_inter in interaction_list]
        elif sort_by == 'RPMAX2':
            rp_list = [max(d_inter._twisted_1, d_inter._twisted_2) \
                       for d_inter in interaction_list]
        else:
            print('[ERROR] Invalid \'sort_by\' parameter! Must be \'RPNUM\' or \'RPMAX\'')

        # Sort interaction indices by read pair number
        sorted_idx = [i for _, i in sorted(zip(rp_list, idx_list), reverse=True)]

        # Select indices of the top n interactions
        sorted_idx = sorted_idx[0:top_n]

        # Create and return new interaction set with top n interactions
        return_interaction_set = DiachromaticInteractionSet()
        for idx in sorted_idx:
            return_interaction_set._inter_dict[interaction_list[idx].key] = interaction_list[idx]
        return return_interaction_set
