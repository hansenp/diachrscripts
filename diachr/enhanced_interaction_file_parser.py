import gzip
import os
from typing import Tuple

from .chromosomal_position import ChromosomalPosition


class EnhancedInteractionFileParser:
    def __init__(self, path):
        if not os.path.exists(path):
            raise ValueError("Could not find file at %s" % path)
        self.path = path
        self.dir_inter_num = 0
        self.undir_inter_num = 0
        self.dir_inter_aa_num = 0
        self.undir_inter_aa_num = 0
        self.dir_inter_aa_rp_array = []

    def _parse_line(self, line):
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
        field = line.rstrip('\n').split("\t")
        # Split string for digest pair associated with the interaction
        coordinate_pair = field[0].split(";")
        chrompos_a = ChromosomalPosition(coordinate_pair[0])
        chrompos_b = ChromosomalPosition(coordinate_pair[1])

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
        rp_simple = int(field[4].split(":")[0]) 
        rp_twisted = int(field[4].split(":")[1])

        i_dist = int(field[1])

        return chrompos_a, chrompos_b, syms_a, tsss_a, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_simple, rp_twisted, i_dist

  
    def parse(self, callbacks: Tuple["ParseCallback"]):
            print("\t[INFO] Iterating enhanced interaction file ...")
            with gzip.open(self.path, 'rt') as fp:
                n_interaction_total = 0
                for line in fp:
                    n_interaction_total += 1
                    if n_interaction_total % 1000000 == 0:
                        print("\t\t[INFO]", n_interaction_total, "interactions processed ...")
                    # Parse enhanced interactions line
                    chrompos_a, chrompos_b, syms_a, tsss_a, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_simple, rp_twisted, i_dist = \
                        self._parse_line(line)
                    for callback in callbacks:
                        callback.on_end_line(chrompos_a, chrompos_b, syms_a, tsss_a, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_simple, rp_twisted, i_dist)
                    #print("chr_a {}, syms_a {} tsss_a {}".format(chrompos_a, syms_a, tsss_a))
                    #print("chr_b {}, syms_b {} tsss_b {}".format(chrompos_b, syms_b, tsss_b))
                    #print("enrichment_pair_tag {} strand_pair_tag {} interaction_category {}".format(enrichment_pair_tag, strand_pair_tag, interaction_category))
                    #print("neg_log_p_value {}, rp_simple {}, rp_twisted {}, i_dist {}".format(neg_log_p_value, rp_simple, rp_twisted, i_dist ))