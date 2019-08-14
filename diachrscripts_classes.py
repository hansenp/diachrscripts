from scipy.stats import binom
import gzip

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


class Interaction:

    # Class to represent a genomic interaction between two restriction digests

    # Class attributes
    digest_distance = None  # Distance between the two interacting digests
    n_simple = 0            # Number of simple read pairs
    n_twisted = 0           # Number of twisted read pairs
    cis = None              # Intrachromosomal interaction
    type = None             # Either simple, twisted or undirected

    # Initializer / instance Attributes
    def __init__(self, digest_1, digest_2, n_simple, n_twisted):
        self.digest_1 = digest_1
        self.digest_2 = digest_2
        if digest_1.get_chromosome() == digest_2.get_chromosome():
            self.cis = True
        else:
            self.cis = False
        self.n_simple = n_simple
        self.n_twisted = n_twisted

    # Instance method
    def get_digest_distance(self):
        if self.digest_distance == None:
            self.digest_distance = self.digest_2.get_start() - self.digest_1.get_end() # distance between digest ends
        return self.digest_distance

    def is_cis(self):
        return self.cis

    def get_status_pair_flag(self):
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

    def get_binomial_p_value(self): # check this function for small simple and twisted read counts
        if self.n_simple < self.n_twisted:
            return 1 - binom.cdf(self.n_twisted-1, self.n_simple + self.n_twisted, 0.5)
        else:
            return 1 - binom.cdf(self.n_simple-1, self.n_simple + self.n_twisted, 0.5)

    def set_interaction_type(self, type):
        if(type != 'S' and type != 'T' and type != 'U'):
            raise Exception('[FATAL] Invalid interaction type. Should be either \'S\', \'T\' or \'U\' but was {}.', type)
        else:
            self.type = type

    def get_interaction_type(self):
        return self.type

    def get_first_digest(self):
        return self.digest_1

    def get_second_digest(self):
        return self.digest_2


class TSSInfo:

    # Class that represents a TSS/gene

    # Attributes
    strand = None
    gene_id = None
    gene_symbol = None
    fpkm = None

    # Initializer
    def __init__(self, strand, gene_id, gene_symbol):
        self.strand = strand
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol

    # Methods
    def set_fpkm(self, fpkm):
        self.fpkm = fpkm


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
        self.tss_info_dict[gene_symbol] = TSSInfo(strand, gene_id, gene_symbol) # use gene_id as key

    # Methods
    def append_TSSInfo(self, strand, gene_id, gene_symbol):
        self.tss_info_dict[gene_symbol] = TSSInfo(strand, gene_id, gene_symbol) # use gene_id as key

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
            print("\tgene_id:" + info.gene_id + "|" + "gene_symbol:" + info.gene_symbol + "|" + "strand:" + info.strand)






class TSSCoordinateMap:

    # Class that represents all coordinates that have one or more TSS

    # Attributes
    tss_coord_dict = None      # dictionary with <chromosome:pos> keys and TSS as values
    annotation_format = None
    annotation_file_path = None

    # Initializer
    def __init__(self, annotation_file_path, format):
        self.annotation_format = format
        self.annotation_file_path = annotation_file_path
        if format == "refGene":
            self.tss_coord_dict = {}
            self.parse_ref_gene_file()

    # Methods
    def parse_ref_gene_file(self):

        print("[INFO] Parsing " + self.annotation_format + " file: " + self.annotation_file_path + " ...")

        # open file
        with gzip.open(self.annotation_file_path) as fp:

            # get first line
            line = fp.readline()

            # iterate file
            while line:

                # parse line
                values = line.split("\t")
                gene_id = values[1]
                chromosome = values[2]
                strand = values[3]
                if strand == "+":
                    position = values[4]
                else:
                    position = values[5]
                gene_symbol = values[12]

                # construct key
                chr_pos_key = chromosome + ":" + position

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
        print("[INFO] Found " + str(n_coord_total) + " coordinates that host TSS for " + str(n_genes_total) + " genes.")
        print("[INFO] Found " + str(n_coord_multiple_genes) + " coordinates that host TSS for more than one gene.")
        print("[INFO] Found " + str(n_coord_multiple_genes_different_strands) + " coordinates that host TSS on different strands.")


    def get_num_of_coords_with_miltiple_genes(self):
        num = 0
        for coord in self.tss_coord_dict.values():
            if 1 < coord.get_num_of_genes():
                num += 1
        return num

    def add_fpkm_values(self, cuffdiff_genes_fpkm_tracking_file):

        # read FPKM values to hash with gene IDs as keys

        fpkm_hash = {}

        with open(cuffdiff_genes_fpkm_tracking_file) as fp:

            line = fp.readline()

            while line:
                values = line.split("\t")
                gene_id = values[3]
                if gene_id == "gene_id":  # skip first line
                    line = fp.readline()
                    continue
                else:
                    if gene_id in fpkm_hash:
                        print "Multiple FPKM for a gene!"
                    else:
                        fpkm_hash[gene_id] = float(values[9])

                line = fp.readline()

        fp.close()

        # Add FPKM to TSSMap
        cnt_no_fpkm = []
        cnt_yes_fpkm = 0
        for key, tss in self.__chr_pos_to_tss.iteritems():
            if tss.get_gene_id() in fpkm_hash:
                tss.set_fpkm(fpkm_hash[tss.get_gene_id()])
                cnt_yes_fpkm += 1
            else:
                cnt_no_fpkm.append(tss.get_gene_id())

        print("[WARNING] No FPKM for " + str(len(cnt_no_fpkm)) + " TSS out of " + str(len(self.__chr_pos_to_tss)) + " in total. List of gene IDs: " + str(cnt_no_fpkm))
