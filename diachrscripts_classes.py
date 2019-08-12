from scipy.stats import binom
import gzip

class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest.

    # Class Attribute
    __chromosome = None
    __start = None
    __end = None
    __active = False

    # Initializer / Instance Attributes
    def __init__(self, chromosome, start, end):
        self.__chromosome = chromosome
        self.__start = start
        self.__end = end

    def set_active(self):
        self.__active = True

    def is_active(self):
        return self.__active

    def get_chromosome(self):
        return self.__chromosome

    def get_start(self):
        return self.__start

    def get_end(self):
        return self.__end


class Interaction:

    # Class to represent a genomic interaction between two restriction digests.

    # Class attributes
    __digest_distance = None  # Distance between the two interacting digests
    __n_simple = 0            # Number of simple read pairs
    __n_twisted = 0           # Number of twisted read pairs
    __cis = None              # Intrachromosomal interaction
    __type = None             # Either simple, twisted or undirected

    # Initializer / instance Attributes
    def __init__(self, digest_1, digest_2, n_simple, n_twisted):
        self.digest_1 = digest_1
        self.digest_2 = digest_2
        if digest_1.get_chromosome() == digest_2.get_chromosome():
            self.__cis = True
        else:
            self.__cis = False
        self.__n_simple = n_simple
        self.__n_twisted = n_twisted

    # Instance method
    def get_digest_distance(self):
        if self.__digest_distance == None:
            self.__digest_distance = self.digest_2.get_start() - self.digest_1.get_end() # distance between digest ends
        return self.__digest_distance

    def is_cis(self):
        return self.__cis

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
        if self.__n_simple < self.__n_twisted:
            return 1 - binom.cdf(self.__n_twisted-1, self.__n_simple + self.__n_twisted, 0.5)
        else:
            return 1 - binom.cdf(self.__n_simple-1, self.__n_simple + self.__n_twisted, 0.5)

    def set_interaction_type(self, type):
        if(type != 'S' and type != 'T' and type != 'U'):
            raise Exception('[FATAL] Invalid interaction type. Should be either \'S\', \'T\' or \'U\' but was {}.', type)
        else:
            self.__type = type

    def get_interaction_type(self):
        return self.__type

    def get_first_digest(self):
        return self.digest_1

    def get_second_digest(self):
        return self.digest_2


class TSS:

    # Attributes
    __chromosome = None
    __position = None
    __strand = None
    __gene_id = None
    __gene_symbol = None
    __fpkm = None

    # Initializer
    def __init__(self, chromosome, position, strand, gene_id, gene_symbol):
        self.__chromosome = chromosome
        self.__position = position
        self.__strand = strand
        self.__gene_id = gene_id
        self.__gene_symbol = gene_symbol

    def get_gene_id(self):
        return self.__gene_id

    def set_fpkm(self, fpkm):
        self.__fpkm = fpkm

    def get_chromosome(self):
        return self.__chromosome

    def get_position(self):
        return self.__position

    def get_strand(self):
        return self.__strand

    def get_gene_symbol(self):
        return self.__gene_symbol

class TSSMap:

    # Attributes
    __chr_pos_to_tss = {}      # Hash with <chromosome:pos> keys and TSS as values
    __format = None            # Until now only the refGene format

    # Initializer
    def __init__(self, annotation_file, format):
        self.__format = format
        if format == "refGene":
            __tss_to_gene_id = self.__parse_ref_gene_file(annotation_file)

    # Parse UCSC refGene file
    def __parse_ref_gene_file(self, annotation_file):

        print("[INFO] Parsing refGene file...")

        cnt_multiple_tss_at_same_position = 0
        cnt_multiple_gene_symols_at_same_position = 0
        with gzip.open(annotation_file) as fp:
            line = fp.readline()

            while line:
                values = line.split("\t")

                gene_id = values[1]
                chromosome = values[2]
                strand = values[3]
                if strand == "+":
                    position = values[4]
                else:
                    position = values[5]

                chr_pos_key = chromosome + ":" + position

                gene_symbol = values[12]

                if chr_pos_key in self.__chr_pos_to_tss:
                    cnt_multiple_tss_at_same_position += 1;
                    #if strand != self.__chr_pos_to_tss[chr_pos_key].get_strand():
                    #if gene_symbol != self.__chr_pos_to_tss[chr_pos_key].get_gene_symbol():
                    if gene_id != self.__chr_pos_to_tss[chr_pos_key].get_gene_id():
                       # print("---------")
                       # print("Gene ID: " + gene_id + "\t" + self.__chr_pos_to_tss[chr_pos_key].get_gene_id())
                       # print("Chromosome: " + chromosome + "\t" + self.__chr_pos_to_tss[chr_pos_key].get_chromosome())
                       # print("Position: " + position + "\t" + self.__chr_pos_to_tss[chr_pos_key].get_position())
                        #print("Strand: " + strand + "\t" + self.__chr_pos_to_tss[chr_pos_key].get_strand())
                        #print("Gene symbol: " + gene_symbol + "\t" + self.__chr_pos_to_tss[chr_pos_key].get_gene_symbol())
                        cnt_multiple_gene_symols_at_same_position += 1;


                tss = TSS(chromosome, position, strand, gene_id, gene_symbol)
                self.__chr_pos_to_tss[chr_pos_key] = tss

                line = fp.readline()

        fp.close()
        print(cnt_multiple_tss_at_same_position)
        print(cnt_multiple_gene_symols_at_same_position)

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
