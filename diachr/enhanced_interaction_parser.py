import gzip
import os
from typing import Tuple, List
from .enhanced_interaction import EnhancedInteraction


class EnhancedInteractionParser:
    """
    This class coordinates the parsing of enhancer interaction file-format files.
    It creates a list of EnhancedInteraction objects, each representing one data line.
    """

    def __init__(self, path:str):
        if not os.path.exists(path):
            raise ValueError("Could not find file at %s" % path)
        self.path = path
        self.dir_inter_num = 0
        self.undir_inter_num = 0
        self.dir_inter_aa_num = 0
        self.undir_inter_aa_num = 0
        self.dir_inter_aa_rp_array = []

    def _split_chromosome(self, chromstring:str) -> Tuple:
        """
        Expect to get a string like this 'chr6:89716777-89716882
        :param str: Expect to get a string like this 'chr6:89716777-89716882
        :return: return a tuple (chr6, 89716777, 89716882)
        """
        pair = chromstring.split(":")
        if len(pair) != 2:
            raise ValueError("Expected two fields but got %d fields from %s" % (len(pair), str))
        chrom = pair[0]
        positions = pair[1].split("-")
        if len(positions) != 2:
            raise ValueError("Expected  two fields but got %d fields from %s" % (len(positions), str))
        start = int(positions[0])
        end = int(positions[1])
        return chrom, start, end

    def _split_chromosomes(self, chrompair:str) -> Tuple:
        """
        Expect to get a string like this 'chr6:89716777-89716882;chr6:89916736-89920828
        :param str: A strting line this 'chr6:89716777-89716882;chr6:89916736-89920828
        :return: a tuple, (chr6, 89716777, 89716882, chr6, 89916736, 89920828)
        """
        pair = chrompair.split(";")
        if len(pair) != 2:
            raise ValueError("Expected two fields but got %d fields from %s" % (len(pair), str))
        chromA, startA, endA = self._split_chromosome(pair[0])
        chromB, startB, endB = self._split_chromosome(pair[1])
        return chromA, startA, endA, chromB, startB, endB

    def _split_pair_on_semicolon(self, strng:str) -> Tuple:
        """
        Expect to get a pair like this: CBFA2T3,LOC100129697;CBFA2T3
        :param str: CBFA2T3,LOC100129697;CBFA2T3
        :return:a tuple ("CBFA2T3,LOC100129697", "CBFA2T3")
        """
        if strng is None or len(strng) == 0 or strng == ";":
            return "", ""
        pair = strng.split(";")
        if len(pair) != 2:
            raise ValueError("Expected two fields but got %d fields from %s" % (len(pair), str))
        return pair[0], pair[1]

    def parse_line(self, line: str) -> EnhancedInteraction:
        # Parse enhanced interaction line
        field = line.rstrip().split("\t")
        # field[0] is like this: 'chr6:89716777-89716882;chr6:89916736-89920828
        chr_a, sta_a, end_a, chr_b, sta_b, end_b = self._split_chromosomes(field[0])
        # Split string for pair of comma separated lists of gene symbols
        syms_a, syms_b = self._split_pair_on_semicolon(field[3])
        # Split string for pair of comma separated lists of TSS
        tsss_a, tsss_b = self._split_pair_on_semicolon(field[8])
        # Extract letter pair tag for enrichment status of digests
        enrichment_pair_tag = field[5]
        # Extract symbol pair tag for TSS orientations on associated digests
        strand_pair_tag = field[7]
        # Extract category of interaction with respect to directionality
        interaction_category = field[2]
        # Get negative decadic logarithm of directionality P-value
        neg_log_p_value = float(field[6])
        read_pair_left = int(field[4].split(":")[0])
        read_pair_right = int(field[4].split(":")[1])
        i_dist = int(field[1])        

        enh_int = EnhancedInteraction(chrA=chr_a, staA=sta_a, endA=end_a, symsA=syms_a,tsssA=tsss_a,
                                      chrB=chr_b, staB=sta_b, endB=end_b, symsB=syms_b, tsssB=tsss_b,
                                      enr_pair_tag=enrichment_pair_tag, strand_pair_tag=strand_pair_tag,
                                        intera_cat=interaction_category, neg_log_p=neg_log_p_value,
                                        read_pair_left=read_pair_left, read_pair_right=read_pair_right,
                                        i_dist=i_dist)
        return enh_int                               

    def parse(self) -> List:
        ei_list = []
        if self.path.endswith(".gz"):
            with gzip.open(self.path, 'rt') as fp:
                n_progress = 0
                for line in fp:
                    n_progress += 1
                    if n_progress % 1000000 == 0:
                        print("\tProcessed " + str(n_progress) + " interactions ...")
                    enh_int = self.parse_line(line)
                    ei_list.append(enh_int)
            return ei_list
        else:
            with open(self.path) as fp:
                n_progress = 0
                for line in fp:
                    n_progress += 1
                    if n_progress % 1000000 == 0:
                        print("\tProcessed " + str(n_progress) + " interactions ...")
                    enh_int = EnhancedInteractionParser.parse_line(line)
                    ei_list.append(enh_int)
            return ei_list