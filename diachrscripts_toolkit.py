import gzip

import numpy
from scipy.stats import binom
import numpy as np
import warnings
warnings.filterwarnings("error")
from decimal import Decimal
from math import log
import sys


#######################################################################################################################

class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest

    # Attributes
    chromosome = None
    start = None
    end = None
    active = False

    # Initializer

    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    # Methods
    def set_active(self):
        self.active = True

    def is_active(self):
        return self.active

    def get_chromosome(self):
        return self.chromosome

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

#######################################################################################################################

class Interaction:

    # Class to represent a genomic interaction between two restriction digests

    # Attributes

    digest_distance = None  # Distance between the two interacting digests
    n_simple = 0            # Number of simple read pairs
    n_twisted = 0           # Number of twisted read pairs
    cis = None              # Intra-chromosomal interaction
    type = None             # Either simple, twisted or undirected
    digest_1 = None
    digest_2 = None
    p_value = None
    digest_pair_flag_original_order = None

    def __init__(self, diachromatic_interaction_line):

        fields = diachromatic_interaction_line.split("\t")
        self.digest_1 = Digest(fields[0], int(fields[1]), int(fields[2]))
        if fields[3] == 'A':
            self.digest_1.set_active()
        self.digest_2 = Digest(fields[4], int(fields[5]), int(fields[6]))
        if fields[7] == 'A':
            self.digest_2.set_active()

        if ':' in fields[8]:  # line regular Diachromatic file
            fields2 = fields[8].split(":")
            self.n_simple = int(fields2[0])
            self.n_twisted = int(fields2[1])
            self.type = "TBD"  # to be determined using binomial P-value
        else: # line from LRT script
            self.n_simple = fields[8]
            self.n_twisted = fields[9]
            self.type = fields[13].rstrip() # has been determined already using LRT
        #self.set_interaction_type(i_type)

        if self.digest_1.get_chromosome() == self.digest_2.get_chromosome():
            self.cis = True
        else:
            self.cis = False

        self.digest_pair_flag_original_order = fields[3] + fields[7]

    # Methods

    def get_digest_pair_flag_original_order(self):
        return self.digest_pair_flag_original_order

    def get_digest_distance(self):
        if not self.is_cis():
            raise Exception("[FATAL] Didest distance is not defined for trans interactions.")
        if self.digest_distance == None:
            self.digest_distance = self.digest_2.get_start() - self.digest_1.get_end()
        return self.digest_distance

    def is_cis(self):
        return self.cis

    def get_digest_status_pair_flag(self):
        if self.digest_1.is_active():
            state_1 = 'A'
        else:
            state_1 = 'I'
        if self.digest_2.is_active():
            state_2 = 'A'
        else:
            state_2 = 'I'
        category = state_1 + state_2
        return sorted(category)[0]+sorted(category)[1]

    def get_binomial_p_value(self):
        if self.p_value != None:
            return self.p_value
        else:
            if self.n_simple < self.n_twisted:
                self.p_value = 1 - binom.cdf(self.n_twisted-1, self.n_simple + self.n_twisted, 0.5)
                return self.p_value
            else:
                self.p_value = 1 - binom.cdf(self.n_simple-1, self.n_simple + self.n_twisted, 0.5)
                return self.p_value

    def set_logsf_binomial_p_value(self): # natural logarithm
        numpy.seterr(divide='ignore')
        if self.p_value != None:
            return self.p_value
        else:
            if self.n_simple < self.n_twisted:
                self.p_value = binom.logsf(self.n_twisted - 1, self.n_simple + self.n_twisted, 0.5)
                return self.p_value
            else:
                self.p_value = binom.logsf(self.n_simple - 1, self.n_simple + self.n_twisted, 0.5)
                return self.p_value


    def set_interaction_type(self, type, p_thresh, n_indefinable_cutoff):
        """
        This function determines whether this interaction is significantly directed, i.e. simple or twisted, using a
        binomial distribution with p = 0.5 and returns a 'S' for simple and 'T' for twisted.
        Otherwise, the function returns an 'U' for undirected or a 'NA', if the sum of twisted and simple is smaller
        than five because interactions with less read pairs cannot be significant using the binomial distribution.
        If the interaction input file was generated by the LRT script, the type of this interaction is set to the
        given type. Note that interactions from LRT cannot have the type 'NA'.

        :return: 'S' (simple),'T' (twisted), 'U' (undirected), 'NA' (not available)
        """
        if(type != 'S' and type != 'T' and type != 'U' and type != "TBD"):
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'TBD' but was " + type + ".")
        else:
            if type == "TBD":
                if self.n_twisted + self.n_simple < n_indefinable_cutoff:
                    type = 'NA' # an interaction with less than 5 read pairs cannot be significant with a binomial P-value threshold of 0.05
                elif self.get_binomial_p_value() <= p_thresh:
                    if self.n_twisted < self.n_simple:
                        type = 'S'
                    else:
                        type = 'T'
                else:
                    type = 'U'
            self.type = type

    def get_interaction_type(self):
        return self.type

    def get_first_digest(self):
        return self.digest_1

    def get_second_digest(self):
        return self.digest_2

    def is_cis_long_range(self, d_dist):
        if not self.is_cis():
            return False
        if self.get_digest_distance() < d_dist:
            return False
        return True

    def get_coord_string(self):
        d1_coords = self.digest_1.get_chromosome() + ":" + str(self.digest_1.get_start()) + "-" + str(self.digest_1.get_end())
        d2_coords = self.digest_2.get_chromosome() + ":" + str(self.digest_2.get_start()) + "-" + str(self.digest_2.get_end())
        return d1_coords + ";" + d2_coords

#######################################################################################################################

class TSSInfo:

    # Class that represents additional information about TSS beyond coordinates

    # Attributes
    strand = None
    gene_id = None
    gene_symbol = None
    fpkm = None
    expression_category = None

    # Initializer

    def __init__(self, strand, gene_id, gene_symbol):
        self.strand = strand
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol

    # Methods

    def set_fpkm(self, fpkm):
        self.fpkm = fpkm

    def set_expression_category(self, category):
        self.expression_category = category

#######################################################################################################################

class TSSCoordinate:

    # Class that represents a genomic coordinate that corresponds to a TSS for one or more genes

    # Attributes

    chromosome = None
    position = None
    tss_info_dict = None   # stores information for all TSS at this position

    # Initializer

    def __init__(self, chromosome, position, strand, gene_id, gene_symbol):
        self.chromosome = chromosome
        self.position = position
        self.tss_info_dict = {}
        self.tss_info_dict[gene_id] = TSSInfo(strand, gene_id, gene_symbol)

    # Methods

    def append_TSSInfo(self, strand, gene_id, gene_symbol):
        self.tss_info_dict[gene_id] = TSSInfo(strand, gene_id, gene_symbol)

    def get_num_of_genes(self):
        return len(self.tss_info_dict)

    def has_multiple_genes(self):
        return 1 < self.get_num_of_genes()

    def has_tss_on_both_strands(self):
        current_strand = None
        for info in self.tss_info_dict.values():
            if current_strand == None:
                current_strand = info.strand
            elif current_strand != info.strand:
                return True
        return False

    def print_tss_info(self):
        if self.has_tss_on_both_strands():
            print("[INFO] Found TSS on both strands for " + str(self.get_num_of_genes()) + " genes at " + self.chromosome + ":" + str(self.position))
        else:
            print("[INFO] Found TSS for " + str(self.get_num_of_genes()) + " genes at " + self.chromosome + ":" + str(self.position))
        for info in self.tss_info_dict.values():
            print("\tgene_id:" + info.gene_id + "|" + "gene_symbol:" + info.gene_symbol + "|" + "strand:" + info.strand + "|" + "fpkm:" + str(info.fpkm))

#######################################################################################################################

class TSSCoordinateMap:

    # Class that represents all coordinates that have one or more TSS

    # Attributes

    tss_coord_dict = None
    annotation_format = None
    annotation_file_path = None

    fpkm_n_zero = None
    fpkm_upper_first_q = None
    fpkm_upper_second_q = None
    fpkm_upper_third_q = None

    # Initializer

    def __init__(self, annotation_file_path, format):
        self.annotation_format = format
        self.annotation_file_path = annotation_file_path
        if format == "refGene":
            self.tss_coord_dict = {}
            self.parse_ref_gene_file()

    # Methods

    def parse_ref_gene_file(self):

        print("[INFO] Parsing " + self.annotation_format + " annotation file: " + self.annotation_file_path + " ...")

        # open file
        with gzip.open(self.annotation_file_path, mode='rt') as fp:

            # get first line
            line = fp.readline()

            # iterate file
            while line:

                # parse line
                values = line.split('\t')
                gene_id = values[12] # for refGene we use gene symbols as gene ID because refSeq gene IDs correspond to transcripts
                chromosome = values[2]
                strand = values[3]
                if strand == "+":
                    position = int(values[4]) + 1 # UCSC start coordinates are 0-based-> convert to 1-based coordinates so as they are compatible with GOPHER's digest map
                else:
                    position = int(values[5])
                gene_symbol = values[12]

                # construct key
                chr_pos_key = chromosome + ":" + str(position)

                # check whether this coordinate has been seen already
                if chr_pos_key in self.tss_coord_dict.keys():
                    # append TSS info to existing coordinate
                    self.tss_coord_dict[chr_pos_key].append_TSSInfo(strand, gene_id, gene_symbol)
                else:
                    # create new coordinate
                    self.tss_coord_dict[chr_pos_key] = TSSCoordinate(chromosome, position, strand, gene_id, gene_symbol)

                # get next line
                line = fp.readline()

        # close file
        fp.close()

    def analyze_coordinates_and_print_report(self):

        # declare count variables
        n_coord_total = 0
        n_genes_total = 0
        n_coord_multiple_genes = 0
        n_coord_multiple_genes_different_strands = 0

        # iterate coordinates and count
        for coord in self.tss_coord_dict.values():
            n_coord_total += 1
            n_genes_total += coord.get_num_of_genes()
            if 1 < coord.get_num_of_genes():
                n_coord_multiple_genes += 1
                if coord.has_tss_on_both_strands():
                    n_coord_multiple_genes_different_strands += 1

        # print report
        print("\t[INFO] Found " + str(n_coord_total) + " coordinates that host TSS for " + str(n_genes_total) + " genes.")
        print("\t[INFO] Found " + str(n_coord_multiple_genes) + " coordinates that host TSS for more than one gene.")
        print("\t[INFO] Found " + str(n_coord_multiple_genes_different_strands) + " coordinates that host TSS on different strands.")


    def get_num_of_coords_with_miltiple_genes(self):
        num = 0
        for coord in self.tss_coord_dict.values():
            if 1 < coord.get_num_of_genes():
                num += 1
        return num

    def parse_cuffdiff_genes_fpkm_tracking_file(self, cuffdiff_genes_fpkm_tracking_file):

        print("[INFO] Parsing genes.fpkm_tracking file of cuffdiff: " + cuffdiff_genes_fpkm_tracking_file + " ...")

        # read FPKM values to temporary dictionary

        tmp_fpkm_dict = {}

        with open(cuffdiff_genes_fpkm_tracking_file) as fp:

            line = fp.readline() # skip first line
            line = fp.readline()
            while line:
                line_fields = line.split("\t")
                gene_id = line_fields[3]
                fpkm = float(line_fields[9])
                tmp_fpkm_dict[gene_id] = fpkm
                line = fp.readline()

        # determine number of zero FPKMs and quartiles
        tmp_fpkm_list = list(tmp_fpkm_dict.values())
        self.fpkm_n_zero = tmp_fpkm_list.count(0)
        tmp_fpkm_list = list(filter(lambda a: a != 0, tmp_fpkm_list)) # remove zero FPKMs
        self.fpkm_upper_first_q = round(np.quantile(tmp_fpkm_list, .25),2)
        self.fpkm_upper_second_q = round(np.quantile(tmp_fpkm_list, .50),2)
        self.fpkm_upper_third_q = round(np.quantile(tmp_fpkm_list, .75),2)

        # add FPKM values to TSSInfo
        genes_without_fpkm = []
        for coord in self.tss_coord_dict.values():
            for tss_info in coord.tss_info_dict.values():
                if tss_info.gene_id in tmp_fpkm_dict:
                    tss_info.set_fpkm(tmp_fpkm_dict[tss_info.gene_id])
                else:
                    genes_without_fpkm.append(tss_info.gene_id)

        if 0 < len(genes_without_fpkm):
            print("\t[WARNING] No FPKM for the following " + str(len(genes_without_fpkm)) + " genes: " + " ".join(genes_without_fpkm))

        # print info about zero FPKM and quartiles
        print ("\t[INFO] There were", self.fpkm_n_zero, "genes with zero FPKM.")
        print ("\t[INFO] Upper FPKM limit of the first quartile:", self.fpkm_upper_first_q)
        print ("\t[INFO] Upper FPKM limit of the second quartile:", self.fpkm_upper_second_q)
        print ("\t[INFO] Upper FPKM limit of the third quartile:", self.fpkm_upper_third_q)

    def set_expression_categories(self, q_limit):

        print("[INFO] Setting expression level categories...")

        if q_limit == 1:
            thresh = self.fpkm_upper_first_q
            print("\t[INFO] Will use upper limit of the first quartile (" + str(self.fpkm_upper_first_q)  + ") as FPKM threshold for inactive genes.")
        elif q_limit == 2:
            thresh = self.fpkm_upper_second_q
            print("\t[INFO] Will use upper limit of the second quartile (" + str(self.fpkm_upper_second_q) + ") as FPKM threshold for inactive genes.")
        elif q_limit == 3:
            thresh = self.fpkm_upper_third_q
            print("\t[INFO] Will use upper limit of the third quartile (" + str(self.fpkm_upper_third_q) + ") as FPKM threshold for inactive genes.")
        else:
            raise Exception("\t[FATAL] Threshold must be an upper limits of the first three quartiles. Should be either 1, 2 or 3 but was " + str(q_limit) + ".")

        none_fpkm = 0
        inactive_category = 0
        active_category = 0

        for coord in self.tss_coord_dict.values():
            for tss_info in coord.tss_info_dict.values():
                if tss_info.fpkm == None:
                    tss_info.set_expression_category(None)
                    none_fpkm += 1
                elif tss_info.fpkm == 0 or tss_info.fpkm < thresh:
                    tss_info.set_expression_category(0)
                    inactive_category +=1
                else:
                    tss_info.set_expression_category(1)
                    active_category += 1

        print("\t[WARNING] There are " + str(none_fpkm) + " without FPKM. Category was set to 'None'.")
        print("\t[INFO] " + str(inactive_category) + " genes were categorized as inactive.")
        print("\t[INFO] " + str(active_category) + " genes were categorized as active.")

    def has_key(self, key):
        return key in self.tss_coord_dict

    def get_coord_category(self, key):

        coord_expression_categories = [] # note that there may be multiple TSS a given coordinate

        if key in self.tss_coord_dict:
            for tss_info in self.tss_coord_dict[key].tss_info_dict.values():
                coord_expression_categories.append(tss_info.expression_category)
        else:
            return -1 # no TSS annotated at this position

        if 1 in coord_expression_categories:
            return 1 # one active coord is enough to make the digest active
        elif 0 in coord_expression_categories:
            return 0
        else:
            return tss_info.gene_id # no FPKM for annotated TSS

    def get_coord_strand(self, coord_key):

        coord_strands = []  # note that there may be multiple TSS a given coordinate

        if coord_key in self.tss_coord_dict:
            for tss_info in self.tss_coord_dict[coord_key].tss_info_dict.values():
                coord_strands.append(tss_info.strand)
        else:
            return -1

        if '-' in coord_strands and '+' in coord_strands:
            return 'd'
        elif '-' in coord_strands:
            return '-'
        elif '+' in coord_strands:
            return '+'
        else:
            raise Exception("[FATAL] Invalid strand symbol. Should be either '-' or '+' but was " + str(coord_strands) + ".")

    def get_coord_symbols(self, coord_key):
        coord_symbols = []
        if coord_key in self.tss_coord_dict:
            for tss_info in self.tss_coord_dict[coord_key].tss_info_dict.values():
                coord_symbols.append(tss_info.gene_symbol)
        else:
            return -1
        return coord_symbols



#######################################################################################################################

class PairKeyDict:

    states = []

    def __init__(self, states):
        self.states = states
        self.pair_dict = {}
        for i in states:
            for j in states:
                key = i + "/" + j
                self.pair_dict[key] = 0

    def print_pair_key_dict(self):
        for i in self.states:
            for j in self.states:
                key = i + "/" + j
                print(key + "\t" + str(self.pair_dict[key]))


#######################################################################################################################

def get_fraction(value, values):
    """
    :param value: A single value.
    :param values: A list of values that contains 'value'.
    :return: A string representing a formatted number with two decimal places.
    """
    if 0 < sum(values):
        return str(Decimal(float(1.0*value/sum(values))).quantize(Decimal(10) ** -2))
    else:
        return str(0.00)

def get_string_formatted_fraction(numerator, denominator):
    if 0 < denominator or 0 > denominator:
        return str(Decimal(float(1.0 * numerator / denominator)).quantize(Decimal(10) ** -2))
    else:
        return str(0.00)

def get_x_inter_binomial_log10_p_value(k, n, x):
    """
    Calculate a P-value using a binomial distribution for observing k out of n interactions with a total number of x
    read pairs that have either zero simple or zero twisted read pairs, given an equal probability of 0.5 for simple
    and twisted read pairs.

    For instance, a single interaction with zero simple and two twisted read pairs corresponds to sequence of two
    independent Bernoulli trials and in each trial the probability for drawing a twisted read pair is 0.5.
    Therefore, the probability to draw two twisted read pairs in a sequence of two trials is 0.5^2,
    and the probability to draw three twisted read pairs in a sequence of three trials is 0.5^3 and so on,
    i.e. p=0.5^x and since we do not distinguish between zero simple and zero twisted interactions, we have 2*0.5^x.

    This function uses a binomial distribution binom(n,k,p) in order to calculate the P-value.

    :param k: Number of interactions with either zero simple or zero twisted read pairs
    :param n: Total number of interactions with a fixed number of x read pairs
    :param x: Total number of read pairs in n interactions
    :return: Negative decadic logarithm of the binomial P-value and error code
    """
    if n==0: # no interaction, no P-value
        return "NA", 1
    #print("-binom.logsf(" + str(k) + ", " + str(n) + ", 2*(0.5 ** " + str(x) + "))/log(10)")
    try:
        return -binom.logsf(k, n, 2*(0.5 ** x))/log(10), 0
    except RuntimeWarning: # underflow
       return -log(sys.float_info.min * sys.float_info.epsilon)/log(10), 2 # return smallest possible float


def calculate_binomial_p_value(n_simple, n_twisted): # conventional with 1 - CDF
    if n_simple < n_twisted:
        p_value = 1 - binom.cdf(n_twisted - 1, n_simple + n_twisted, 0.5)
        return p_value
    else:
        p_value = 1 - binom.cdf(n_simple - 1, n_simple + n_twisted, 0.5)
        return p_value

def calculate_binomial_logsf_p_value(n_simple, n_twisted): # (natural) logsf
    try:
        if n_simple < n_twisted:
            p_value = binom.logsf(n_twisted - 1, n_simple + n_twisted, 0.5)
            return p_value
        else:
            p_value = binom.logsf(n_simple - 1, n_simple + n_twisted, 0.5)
            return p_value
    except RuntimeWarning:
        return log(sys.float_info.min * sys.float_info.epsilon) # return natural log of smallest possible float

def find_indefinable_n(p_value_cutoff, verbose = True):
    """
    :param p_value_cutoff: Chosen threshold
    :return: Minimal n that yields a significant P-value given the threshold
    """
    if verbose:
        print("[INFO] Looking for smallest number of read pairs n that yields a significant P-value with the given threshold of " + str(p_value_cutoff) + ".")
    for n in range(1, 1000):
        p_value = calculate_binomial_p_value(n, 0)
        if p_value < p_value_cutoff:
            if verbose:
                print("\t[INFO] Smallest n: " + str(n) + " read pairs (" + str(p_value) + ")")
            return n, p_value

def get_n_dict(diachromatic_interaction_file, status_pair_flag, min_digest_dist, p_value_cutoff):
    """
    :param diachromatic_interaction_file: Interaction file in Diachromatic format
    :param status_pair_flag: ALL, AA, AI or II
    :param min_digest_dist: Minimal allowed distance between interacting digests, typically 10,000
    :return: Dictionary with total read pair numbers (n) as keys and the corresponding numbers interactions with n read pairs
    """

    p_val_dict = {} # dictionary that stores previously calculated P-values (saves time)
    n_dict = {}  # dictionary that stores the numbers of interactions with n read pairs
    digest_distances = []

    print("[INFO] Determining distribution of n in significant interactions in: " + diachromatic_interaction_file + " ...")
    with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
        line = fp.readline()
        n_interaction_total = 0
        while line:

            n_interaction_total += 1
            if n_interaction_total % 1000000 == 0:
                print("\t[INFO]", n_interaction_total, "interactions processed ...")

            # parse line representing one interaction
            interaction = Interaction(line)

            # restrict analysis to interactions between targeted promoters
            if status_pair_flag != "ALL" and interaction.get_digest_status_pair_flag() != status_pair_flag:
                line = fp.readline()
                continue

            # restrict analysis to cis long range interactions
            if not (interaction.is_cis_long_range(min_digest_dist)):
                line = fp.readline()
                continue

            n_total = interaction.n_simple + interaction.n_twisted

            key = ":".join((str(interaction.n_simple),str(interaction.n_twisted)))
            if key in p_val_dict:
                p_val = p_val_dict[key]
            else:
                p_val = interaction.get_binomial_p_value()
                p_val_dict[key] = p_val


            if p_val <= p_value_cutoff:
                digest_distances.append(interaction.get_digest_distance())
                if n_total in n_dict:
                    n_dict[n_total] += 1
                else:
                    n_dict[n_total] = 1

            line = fp.readline()

    if(0<len(digest_distances)):
        n_dict[-1] = np.quantile(digest_distances, .25)
        n_dict[-3] = np.quantile(digest_distances, .75)

    print("... done.")

    fp.close()
    return n_dict


def get_n_dict_definable(diachromatic_interaction_file, status_pair_flag, min_digest_dist, n_indef):
    """
    :param diachromatic_interaction_file: Interaction file in Diachromatic format
    :param status_pair_flag: ALL, AA, AI or II
    :param min_digest_dist: Minimal allowed distance between interacting digests, typically 10,000
    :return: Dictionary with total read pair numbers (n) as keys and the corresponding numbers interactions with n read pairs
    """

    p_val_dict = {} # dictionary that stores previously calculated P-values (saves time)
    n_dict = {}  # dictionary that stores the numbers of interactions with n read pairs
    digest_distances = []

    print("[INFO] Determining distribution of n in significant interactions in: " + diachromatic_interaction_file + " ...")
    with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
        line = fp.readline()
        n_interaction_total = 0
        while line:

            n_interaction_total += 1
            if n_interaction_total % 1000000 == 0:
                print("\t[INFO]", n_interaction_total, "interactions processed ...")

            # parse line representing one interaction
            interaction = Interaction(line)

            # restrict analysis to interactions between targeted promoters
            if status_pair_flag != "ALL" and interaction.get_digest_status_pair_flag() != status_pair_flag:
                line = fp.readline()
                continue

            # restrict analysis to cis long range interactions
            if not (interaction.is_cis_long_range(min_digest_dist)):
                line = fp.readline()
                continue

            n_total = interaction.n_simple + interaction.n_twisted

            key = ":".join((str(interaction.n_simple),str(interaction.n_twisted)))
            if key in p_val_dict:
                p_val = p_val_dict[key]
            else:
                p_val = interaction.get_binomial_p_value()
                p_val_dict[key] = p_val


            if n_total >= n_indef:
                digest_distances.append(interaction.get_digest_distance())
                if n_total in n_dict:
                    n_dict[n_total] += 1
                else:
                    n_dict[n_total] = 1

            line = fp.readline()

    if(0<len(digest_distances)):
        n_dict[-1] = np.quantile(digest_distances, .25)
        n_dict[-3] = np.quantile(digest_distances, .75)

    print("... done.")

    fp.close()
    return n_dict


def parse_enhanced_interaction_line_with_gene_symbols(interaction_line_with_gene_symbols):
    """
    This function takes a tab separated line with coordinates and gene symbols that correspond to an interaction and
    parses it into individual fields that are relevant for the analyses performed in this script.

    :param interaction_line_with_gene_symbols: Line from a file that was created using the script 'get_gene_symbols_interactions'.
    :return: chr_a: Chromosome of the first digest of the interaction (e.g. chr1)
             sta_a: Start coordinate of the first digest of the interaction (e.g. 123456)
             end_a: End coordinate of the first digest of the interaction (e.g. 234567)
             syms_a: Comma separated list of gene symbols associated with TSS on the first digest (e.g. HIST1H1E,HIST1H2BD)
             tsss_a: Comma separated list of TSS coordinates on the first digest (e.g. chr6:26158121:+,chr6:26156331:+)
             chr_b, sta_b, end_b, syms_b, tsss_b: Same as for the first digest but for the second digest
             enrichment_pair_tag: Two letter tag indicating the enrichment status of the two digests (e.g.  AA, AI, II)
             strand_pair_tag: Two symbol tag separated by '/' indicating the strands of TSS on first and second digest (e.g. -/-, +/+, -/+, +/-, -/d, ...)
             interaction_category: Interaction category with respect to directionality (S, T, URAA, URAI, ...)
    """

    # Split line into individual fields
    field = interaction_line_with_gene_symbols.split("\t")

    # Split string for digest pair associated with the interaction
    coordinate_pair = field[0].split(";")
    coordinate_a = coordinate_pair[0]
    coordinate_b = coordinate_pair[1]

    # Extract digest coordinates
    chr_a = coordinate_a.split(":")[0]
    chr_b = coordinate_b.split(":")[0]

    sta_a = int(coordinate_a.split(":")[1].split("-")[0])
    sta_b = int(coordinate_b.split(":")[1].split("-")[0])

    end_a = int(coordinate_a.split(":")[1].split("-")[1])
    end_b = int(coordinate_b.split(":")[1].split("-")[1])

    # Split string for pair of comma separated lists of gene symbols
    symbols_pair = field[3].split(";")
    syms_a = symbols_pair[0]
    syms_b = symbols_pair[1]

    # Split string for pair of comma separated lists of TSS
    tsss_pair = field[8].split(";")
    tsss_a = tsss_pair[0]
    tsss_b = tsss_pair[1]

    # Extract letter pair tag for enrichment status of digests
    enrichment_pair_tag = field[5]

    # Extract symbol pair tag for TSS orientations on associated digests
    strand_pair_tag = field[7]

    # Extract category of interaction with respect to directionality
    interaction_category = field[2]

    # Get negative decadic logarithm of directionality P-value
    neg_log_p_value = float(field[6])

    # Get total number of read pairs
    rp_total = int(field[4].split(":")[0]) + int(field[4].split(":")[1])

    i_dist = int(field[1])

    return chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist

def set_interaction_category_in_enhanced_interaction_line(line, new_category):
    fields = line.rstrip("\n").split("\t")
    new_line = fields[0] + "\t" + fields[1] + "\t" + new_category + "\t" + fields[3] + "\t" + fields[4] + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8]
    return new_line