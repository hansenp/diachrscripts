import gzip
import os
import copy
from numpy import exp, log
from collections import defaultdict
from .diachromatic_interaction import DiachromaticInteraction
from .diachromatic_interaction import DiachromaticInteraction11
from diachr.binomial_model import BinomialModel


class DiachromaticInteractionSet:
    """
    This class can read interactions from one ore more Diachromatic interaction files and perform the following
    operations on them:

        1. Calculation of P-values and categorization into directed and undirected
        2. Selection of undirected reference interactions
        3. Writing interactions that occur in a required number of input files to a new file

    If interactions have been evaluated and categorized (1,2), the output format is expanded by two columns on the left.
    Column 10 contains the negative of the natural logarithm of the P-value and column 11 the interaction category
    (DI, UI or UIR).

    The functions of this class are tested and detailed documentation on this class can be found in the relevant
    section in the RTD of this repository.
    """

    def __init__(self, enriched_digests_file: str = None):

        # Dictionary that contains all interaction objects
        self._inter_dict = defaultdict(DiachromaticInteraction)

        # Object for calculating P-values
        self._p_values = BinomialModel()

        # Dictionary with information about read files
        self._read_file_info_dict = {'I_FILE': [], 'I_NUM': []}

        # Dictionary with information about the last writing process
        self._write_file_info_dict = {}

        # Dictionary with information about the last evaluation and categorization process
        self._eval_cat_info_dict = {}

        # Dictionary with information about the last reference selection process
        self._select_ref_info_dict = {}

        # Take information about digests selected for enrichment from BED file
        self._enriched_digests_set = None
        if enriched_digests_file is not None:
            print("[INFO] Reading list with digests selected for enrichment ...")
            self._enriched_digests_set = set()
            with open(enriched_digests_file, 'rt') as fp:
                for line in fp:
                    chr, sta, end = line.rstrip().split('\t')
                    self._enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))
            print("\t[INFO] Read " + str(len(self._enriched_digests_set)) + " digests ...")
            print("[INFO] ... done.")

    def parse_file(self, i_file: str = None, verbose: bool = False):
        """
        Parses a file with interactions. For interactions that have already been parsed from another file,
        only the simple and twisted read pair counts are added.

        :param i_file: File with interactions
        """

        if not os.path.exists(i_file):
            raise ValueError("Could not find file at %s" % i_file)

        if verbose:
            print("[INFO] Parsing Diachromatic interaction file ...")
            print("\t[INFO] " + i_file)

        if i_file.endswith(".gz"):
            with gzip.open(i_file, 'rt') as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if verbose:
                        if n_lines % 1000000 == 0:
                            print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.n_simple,
                                                                              twisted=d_inter.n_twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter
        else:
            with open(i_file) as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if verbose:
                        if n_lines % 1000000 == 0:
                            print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.n_simple,
                                                                              twisted=d_inter.n_twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter

        # Keep track of the name of the input file and the number of interactions
        self._read_file_info_dict['I_FILE'].append(i_file)
        self._read_file_info_dict['I_NUM'].append(n_lines)

        if verbose:
            print("[INFO] ... done.")

    def _parse_line(self, line: str = None) -> DiachromaticInteraction:
        """
        Parses a Diachromatic interaction formatted line with 9 fields.

        :param line: Tab separated string 9 fields.
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
        st_string = F[8]  # something like 2:1, representing S:T, simple and twisted counts
        st_fields = st_string.split(":")
        if len(st_fields) != 2:
            raise ValueError("Malformed simple:twisted field in diachromatic line: " + line)
        simple = int(st_fields[0])
        twisted = int(st_fields[1])

        # Get enrichment status of digests from file for enriched digests if available
        if self._enriched_digests_set != None:
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

        di_inter = DiachromaticInteraction(
            chrA=chrA,
            fromA=fromA,
            toA=toA,
            statusA=statusA,
            chrB=chrB,
            fromB=fromB,
            toB=toB,
            statusB=statusB,
            simple=simple,
            twisted=twisted)

        return di_inter

    def evaluate_and_categorize_interactions(self, nln_pval_thresh: float = None, verbose: bool = False):
        """
        Calculate the P-value and define interaction category ('DI' or 'UI') for all interactions in this object.
        'DiachromaticInteraction' objects will be replaced by 'DiachromaticInteraction11' objects.
        """

        if verbose:
            print("[INFO] Rate and categorize interactions ...")

        # Get smallest number of read pairs required for significance
        min_rp, min_rp_pval = self._p_values.find_smallest_significant_n(exp(-nln_pval_thresh), verbose=False)

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
                    print("\t[INFO] Processed " + str(n_processed) + " interactions ...")

            # Discard interactions that don't have enough read pairs to be significant
            if d_inter.rp_total < min_rp:
                n_discarded += 1
                continue

            # Get P-value
            nln_p_val = self._p_values.get_binomial_nnl_p_value(d_inter._simple, d_inter._twisted)

            # Determine interaction category
            if nln_pval_thresh <= nln_p_val:
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
                simple=d_inter.n_simple,
                twisted=d_inter.n_twisted,
                nln_pval=nln_p_val)

            # Set interaction category
            di_11_inter.set_category(i_category)
            d11_inter_dict[d_inter.key] = di_11_inter

        # Replace old list with new list of interactions
        self._inter_dict = d11_inter_dict

        # Prepare dictionary for report
        report_dict = {
            'NLN_PVAL_THRESH': [nln_pval_thresh],
            'MIN_RP': [min_rp],
            'MIN_RP_PVAL': [min_rp_pval],
            'N_PROCESSED': [n_processed],
            'N_DISCARDED': [n_discarded],
            'N_UNDIRECTED': [n_ui],
            'N_DIRECTED': [n_di]
        }

        if verbose:
            print("[INFO] ...done.")

        self._eval_cat_info_dict = report_dict
        return report_dict

    def select_reference_interactions(self, verbose: bool = False):
        """
        Select reference interactions that match directed interactions in terms of enrichment category and total number
        of read pairs per interaction and return a dictionary with information on this selection process.

        :return: Dictionary with information on this selection process
        """

        if verbose:
            print("[INFO] Select reference interactions ...")

        # Nested dictionary that stores the numbers of interactions (value) for different read pair numbers (key)
        rp_inter_dict = {'NN': {},
                         'NE': {},
                         'EN': {},
                         'EE': {}}

        if verbose:
            print("\t[INFO] First pass: Count directed interactions for different read pair counts ...")
        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() == 'DI':

                # Get enrichment status tag pair and read pair number
                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                rp_total = d11_inter.rp_total

                if rp_total not in rp_inter_dict[enrichment_pair_tag]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] = 1
                else:
                    rp_inter_dict[enrichment_pair_tag][rp_total] += 1

        rp_inter_dict_before = copy.deepcopy(rp_inter_dict)
        ui_inter_dict = {'NN': 0,
                         'NE': 0,
                         'EN': 0,
                         'EE': 0}

        if verbose:
            print("\t[INFO] Second pass: Select undirected reference interactions for different read pair counts ...")
        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() != 'DI':

                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                rp_total = d11_inter.rp_total

                if rp_total in rp_inter_dict[enrichment_pair_tag] and 0 < rp_inter_dict[enrichment_pair_tag][rp_total]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] -= 1
                    d11_inter.set_category('UIR')
                else:
                    ui_inter_dict[enrichment_pair_tag] += 1

        # Prepare dictionary for report
        report_dict = {}
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            report_dict['DI_' + enr_cat] = [sum(rp_inter_dict_before[enr_cat].values())]
            report_dict['UIR_' + enr_cat] = [
                sum(rp_inter_dict_before[enr_cat].values()) - sum(rp_inter_dict[enr_cat].values())]
            report_dict['M_UIR_' + enr_cat] = [sum(rp_inter_dict[enr_cat].values())]
            report_dict['UI_' + enr_cat] = [ui_inter_dict[enr_cat]]

        if verbose:
            print("[INFO] ...done.")

        self._select_ref_info_dict = report_dict
        return report_dict

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
        n_has_required_data = 0
        n_incomplete_data = 0
        for d_inter in self._inter_dict.values():
            if d_inter.has_data_for_required_replicate_num(required_replicates):
                n_has_required_data += 1
                out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
            else:
                n_incomplete_data += 1
        out_fh.close()

        # Prepare dictionary for report
        report_dict = {
            'TARGET_FILE': [target_file],
            'INTERACTIONS_NUMBERS': [self._read_file_info_dict['I_NUM']],
            'REQUIRED_REPLICATES': [required_replicates],
            'HAS_ALL_DATA': [n_has_required_data],
            'INCOMPLETE_DATA': [n_incomplete_data]
        }

        if verbose:
            print("[INFO] ... done.")

        self._write_file_info_dict = report_dict
        return report_dict

    def get_read_file_info_dict(self):
        """
        :return: Dictionary that contains information about parsed interaction data.
        """

        read_file_info_dict = copy.deepcopy(self._read_file_info_dict)
        read_file_info_dict['I_FILE'].append('UNION')
        read_file_info_dict['I_NUM'].append(len(self._inter_dict))

        return read_file_info_dict

    def get_read_file_info_report(self):
        """
        :return: String that contains information about parsed interaction data.
        """

        file_num = len(self._read_file_info_dict['I_FILE'])
        report = "[INFO] Report on reading files:" + '\n'
        report += "\t[INFO] Read interaction data from " + str(file_num) + " files:" + '\n'
        for i in range(0, file_num):
            report += "\t\t[INFO] " + str(self._read_file_info_dict['I_NUM'][i]) + " interactions from " + \
                      self._read_file_info_dict['I_FILE'][i] + '\n'
        report += "\t[INFO] The union of all interactions has " + \
                  str(len(self._inter_dict)) + " interactions." + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_write_file_info_dict(self):
        """
        :return: Dictionary that contains information about the last writing process.
        """
        return self._write_file_info_dict

    def get_write_file_info_report(self):
        """
        :return: String that contains information about the last writing process.
        """

        report = "[INFO] Report on writing files:" + '\n'
        report += "\t[INFO] Wrote interactions that occur in at least " + \
                  str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + \
                  " replicates to: " + self._write_file_info_dict['TARGET_FILE'][0] + '\n'
        report += "\t[INFO] Interactions that occur in at least " + \
                  str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + " replicates: " + \
                  str(self._write_file_info_dict['HAS_ALL_DATA'][0]) + '\n'
        report += "\t[INFO] Other interactions: " + \
                  str(self._write_file_info_dict['INCOMPLETE_DATA'][0]) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_write_file_info_table_row(self):
        """
        :return: String consisting of a header line and a line with values relating to written interactions
        """

        table_row = "TARGET_FILE" + "\t" + \
                    "I_NUMS" + "\t" + \
                    "I_NUM_UNION" + "\t" + \
                    "REQUIRED_REPLICATES" + "\t" + \
                    "HAS_ALL_DATA" + "\t" + \
                    "INCOMPLETE_DATA" + '\n'
        table_row += str(self._write_file_info_dict['TARGET_FILE'][0]) + "\t" + \
                     str(self._read_file_info_dict['I_NUM']) + "\t" + \
                     str(len(self._inter_dict)) + "\t" + \
                     str(self._write_file_info_dict['REQUIRED_REPLICATES'][0]) + "\t" + \
                     str(self._write_file_info_dict['HAS_ALL_DATA'][0]) + "\t" + \
                     str(self._write_file_info_dict['INCOMPLETE_DATA'][0]) + '\n'

        return table_row

    def get_eval_cat_info_dict(self):
        """
        :return: Dictionary that contains information about the last evaluation and categorization process.
        """
        return self._eval_cat_info_dict

    def get_eval_cat_info_report(self):
        """
        :return: String that contains information about the last evaluation and categorization process.
        """

        report = "[INFO] Report on evaluation and categorization interactions:" + '\n'
        report += "\t[INFO] Minimum number of read pairs required for significance: " + str(
            self._eval_cat_info_dict['MIN_RP'][0]) + '\n'
        report += "\t[INFO] Corresponding largest P-value: " + "{:.2f}".format(
            -log(self._eval_cat_info_dict['MIN_RP_PVAL'][0])) + '\n'
        report += "\t[INFO] Processed interactions: " + str(self._eval_cat_info_dict['N_PROCESSED'][0]) + '\n'
        report += "\t[INFO] Discarded interactions: " + str(self._eval_cat_info_dict['N_DISCARDED'][0]) + '\n'
        report += "\t[INFO] Not significant interactions (UI): " + str(
            self._eval_cat_info_dict['N_UNDIRECTED'][0]) + '\n'
        report += "\t[INFO] Significant interactions (DI): " + str(self._eval_cat_info_dict['N_DIRECTED'][0]) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_eval_cat_info_table_row(self, out_prefix: str = None):
        """
        :return: String consisting of a header line and a line with values relating to evaluation and categorization
        of interactions
        """

        table_row = "OUT_PREFIX" + '\t' + \
                    "NLN_PVAL_THRESH" + '\t' + \
                    "MIN_RP" + '\t' + \
                    "MIN_RP_NLN_PVAL" + '\t' + \
                    "N_PROCESSED" + '\t' + \
                    "N_DISCARDED" + '\t' + \
                    "N_UNDIRECTED" + '\t' \
                                     "N_DIRECTED" + '\n'

        table_row += str(out_prefix) + '\t' + \
                     "{:.2f}".format(self._eval_cat_info_dict['NLN_PVAL_THRESH'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['MIN_RP'][0]) + '\t' + \
                     "{:.2f}".format(-log(self._eval_cat_info_dict['MIN_RP_PVAL'][0])) + '\t' + \
                     str(self._eval_cat_info_dict['N_PROCESSED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_DISCARDED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_UNDIRECTED'][0]) + '\t' + \
                     str(self._eval_cat_info_dict['N_DIRECTED'][0]) + '\n'

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

        report = "[INFO] Report on selection of undirected reference interactions:" + '\n'
        report += "\t[INFO] Numbers of directed interactions" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict['DI_' + enr_cat][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + str(
                self._select_ref_info_dict['DI_' + enr_cat][0]) + '\n'
        report += "\t\t[INFO] Total: " + str(total) + '\n'
        report += "\t[INFO] Numbers of undirected reference interactions" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict['UIR_' + enr_cat][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + str(
                self._select_ref_info_dict['UIR_' + enr_cat][0]) + '\n'
        report += "\t\t[INFO] Total: " + str(total) + '\n'
        report += "\t[INFO] Numbers of missing undirected reference interactions" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict['M_UIR_' + enr_cat][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + str(
                self._select_ref_info_dict['M_UIR_' + enr_cat][0]) + '\n'
        report += "\t\t[INFO] Total: " + str(total) + '\n'
        report += "\t[INFO] Numbers undirected interactions" + '\n'
        total = 0
        for enr_cat in ['NN', 'NE', 'EN', 'EE']:
            total += self._select_ref_info_dict['UI_' + enr_cat][0]
            report += "\t\t[INFO] Interactions in " + enr_cat + ": " + str(
                self._select_ref_info_dict['UI_' + enr_cat][0]) + '\n'
        report += "\t\t[INFO] Total: " + str(total) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_select_ref_info_table_row(self, out_prefix: str = None):
        """
        :return: String consisting of a header line and a line with values relating to the last selection of reference
        interactions
        """

        table_row = "OUT_PREFIX" + '\t' + \
 \
                    "DI_NN" + '\t' + \
                    "DI_NE" + '\t' + \
                    "DI_EN" + '\t' + \
                    "DI_EE" + '\t' + \
 \
                    "UIR_NN" + '\t' + \
                    "UIR_NE" + '\t' + \
                    "UIR_EN" + '\t' + \
                    "UIR_EE" + '\t' + \
 \
                    "M_UIR_NN" + '\t' + \
                    "M_UIR_NE" + '\t' + \
                    "M_UIR_EN" + '\t' + \
                    "M_UIR_EE" + '\t' + \
 \
                    "UI_NN" + '\t' + \
                    "UI_NE" + '\t' + \
                    "UI_EN" + '\t' + \
                    "UI_EE" + '\n'

        table_row += str(out_prefix) + '\t' + \
 \
                     str(self._select_ref_info_dict["DI_NN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["DI_NE"][0]) + '\t' + \
                     str(self._select_ref_info_dict["DI_EN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["DI_EE"][0]) + '\t' + \
 \
                     str(self._select_ref_info_dict["UIR_NN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UIR_NE"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UIR_EN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UIR_EE"][0]) + '\t' + \
 \
                     str(self._select_ref_info_dict["M_UIR_NN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["M_UIR_NE"][0]) + '\t' + \
                     str(self._select_ref_info_dict["M_UIR_EN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["M_UIR_EE"][0]) + '\t' + \
 \
                     str(self._select_ref_info_dict["UI_NN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UI_NE"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UI_EN"][0]) + '\t' + \
                     str(self._select_ref_info_dict["UI_EE"][0]) + '\n'

        return table_row

    @property
    def interaction_list(self):
        return self._inter_dict.values()
