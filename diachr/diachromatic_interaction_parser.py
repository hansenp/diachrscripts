import gzip
import os
from collections import defaultdict
from .diachromatic_interaction import DiachromaticInteraction
from .diachromatic_interaction import DiachromaticInteraction11
from diachr.binomial_model import BinomialModel

class DiachromaticInteractionParser:
    """
    This class coordinates the parsing of Diachromatic interaction files.
    It creates a list of DiachromaticInteraction objects, each representing one data line.
    """

    def __init__(self, enriched_digests_file: str = None):

        # Dictionary for read files
        self._file_dict = {}

        # Dictionary that contains all interaction objects
        self._inter_dict = defaultdict(DiachromaticInteraction)

        # Dictionary that contains all interaction objects
        self._p_values = BinomialModel()

        # Take information about digests selected for enrichment from BED file
        self._enriched_digests_set = None
        if enriched_digests_file != None:
            print("[INFO] Reading list with digests selected for enrichment ...")
            self._enriched_digests_set = set()
            with open(enriched_digests_file, 'rt') as fp:
                for line in fp:
                    chr, sta, end = line.rstrip().split('\t')
                    self._enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))
            print("\t[INFO] Read " + str(len(self._enriched_digests_set)) + " digests ...")
            print("[INFO] ... done.")


    def parse_file(self, i_file):

        """
        Parses a file with interactions. For interactions that have already been parsed from another file,
        only the simple and twisted read pair counts are added.
        :param i_file: File with interactions
        """

        if not os.path.exists(i_file):
            raise ValueError("Could not find file at %s" % i_file)

        print("[INFO] Parsing Diachromatic interaction file ...")

        print("\t[INFO] " + i_file)

        if i_file.endswith(".gz"):
            with gzip.open(i_file, 'rt') as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if n_lines % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.n_simple, twisted=d_inter.n_twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter
        else:
            with open(i_file) as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if n_lines % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.n_simple, twisted=d_inter.n_twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter

        self._file_dict[i_file] = n_lines

        print("[INFO] ... done.")


    def _parse_line(self, line: str) -> DiachromaticInteraction:
        """
        Parses a Diachromatic interaction formatted line with 9 fields.
        :param line: Diachromatic interaction formatted line with 9 fields.
        :return: DiachromaticInteraction object
        """

        F = line.rstrip().split('\t')
        if len(F) < 9:
            raise ValueError("Malformed line {}".format(line))

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


    def get_file_dict_info(self):
        """
        :return: String that contains information about parsed data.
        """
        file_dict_info = "[INFO] Read interaction data from " + str(len(self._file_dict)) + " files:" + '\n'
        for k, v in self._file_dict.items():
            file_dict_info += "\t[INFO] " + str(v) + " interactions from " + k  + '\n'
        file_dict_info += "[INFO] The union of all interactions has "  + str(len(self._inter_dict)) + " interactions." + '\n'
        return(file_dict_info)

    def write_diachromatic_interaction_file(self, required_replicates = 1, target_file_name = None):
        """
        Writes interactions that occur in a specified minimum number of replicates to a file
        and returns a string with information on this writing process.

        :param required_replicates: Minimum number of replicates
        :param target_file_name: Generated file with interactions
        :return: String with information on this writing process
        """
        write_info = "[INFO] Writing interactions that occur in at least " + str(required_replicates) + " replicates to: " + target_file_name + '\n'
        out_fh = gzip.open(target_file_name, 'wt')
        n_has_required_data = 0
        n_incomplete_data = 0
        for d_inter in self._inter_dict.values():
            if d_inter.has_data_for_required_replicate_num(required_replicates):
                n_has_required_data += 1
                out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
            else:
                n_incomplete_data += 1
        out_fh.close()
        write_info += "\t[INFO] Interactions that occur in at least " + str(required_replicates) + " replicates: " + str(n_has_required_data) + '\n'
        write_info += "\t[INFO] Other interactions: " + str(n_incomplete_data) + '\n'
        write_info += '\n'
        write_info += "FILE_NAME" + "\t" + "INTERACTIONS_NUMBERS" + "\t" + "REQUIRED_REPLICATES" + "\t" + "HAS_ALL_DATA" + "\t" + "INCOMPLETE_DATA" + '\n'
        write_info += str(target_file_name) + "\t" + str(list(self._file_dict.values())) + "\t" + str(required_replicates) + "\t" + str(n_has_required_data) + "\t" + str(n_incomplete_data) + '\n'
        return write_info

    @property
    def file_num(self):
        return len(self._file_dict)

    @property
    def interaction_list(self):
        return self._inter_dict.values()


    def rate_and_categorize_interactions(self, nln_pval_thresh: float):
        """
        Calculates the P-value and defines categories for all interactions in this object.
        DiachromaticInteraction objects will be replaced by DiachromaticInteraction11 objects.
        """

        print("[INFO] Rate and categorize interactions ...")

        d11_inter_dict = {}
        for d_inter in self._inter_dict.values():

            # Get P-value
            nln_p_val = self._p_values.get_binomial_nnl_p_value(d_inter._simple, d_inter._twisted)

            # Determine interaction category
            if nln_pval_thresh <= nln_p_val:
                i_category = "DI"
            else:
                i_category = "UI"

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

            di_11_inter.set_category(i_category)

            d11_inter_dict[d_inter.key] = di_11_inter
            self._inter_dict = d11_inter_dict

        print("[INFO] ...done.")


    def select_reference_interactions(self):
        """
        XXX
        """

        print("[INFO] Select reference interactions ...")

        # Nested dictionary that stores the numbers of interactions (value) for different read pair numbers (key)
        rp_inter_dict = {'NN': {},
                         'NE': {},
                         'EN': {},
                         'EE': {}}

        # First pass: Count directed interactions for different read pair counts
        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() == 'DI':

                # Get enrichment status tag pair and read pair number
                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                rp_total = d11_inter.rp_total

                if rp_total not in rp_inter_dict[enrichment_pair_tag]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] = 1
                else:
                    rp_inter_dict[enrichment_pair_tag][rp_total] += 1


        # Second pass: Select undirected reference interactions for different read pair counts
        for d11_inter in self._inter_dict.values():

            if d11_inter.get_category() != 'DI':

                enrichment_pair_tag = d11_inter.enrichment_status_tag_pair
                rp_total = d11_inter.rp_total

                if rp_total in rp_inter_dict[enrichment_pair_tag] and 0 < rp_inter_dict[enrichment_pair_tag][rp_total]:
                    rp_inter_dict[enrichment_pair_tag][rp_total] -= 1
                    d11_inter.set_category('UIR')

        # Check whether there are read pair numbers for which no undirected reference interactions were found
        for enr_cat in rp_inter_dict.keys():
            print(enr_cat)
            n_missing = 0
            for key in rp_inter_dict[enr_cat].keys():
                if rp_inter_dict[enr_cat][key] < 0:
                    n_missing += rp_inter_dict[enr_cat][key]
            print("Missing interactions: " + str(n_missing))

        print("[INFO] ...done.")
