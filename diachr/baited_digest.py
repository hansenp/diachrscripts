from .diachromatic_interaction import DiachromaticInteraction
from .diachromatic_interaction import DiachromaticInteraction11
import numpy as np


class BaitedDigest:
    """
    Objects of this class are used to group interactions interactions that end in the same bait. Interactions are broken
    down by category and direction as seen from the bait. For the interaction categories directed (DI), undirected
    reference (UIR), undirected (UI) and all interactions (ALL) there are two lists of interactions each, one with
    interactions from the bait to the left (NE) and one with interactions from the bait to the right (EN).
    """

    # Dictionary with lists of interactions ending in baited digest
    interactions = None

    curb_nums = None

    def __init__(self):

        self.interactions = {
            'DI': {
                'NE': [],
                'EN': []},
            'UIR': {
                'NE': [],
                'EN': []},
            'UI': {
                'NE': [],
                'EN': []},
            'ALL': {
                'NE': [],
                'EN': []}
        }

        self.curb_nums = {
            'DI': {
                'NE': 0,
                'EN': 0},
            'UIR': {
                'NE': 0,
                'EN': 0},
            'UI': {
                'NE': 0,
                'EN': 0},
            'ALL': {
                'NE': 0,
                'EN': 0}
        }

    def add_interaction(self, d11_inter: DiachromaticInteraction11):
        self.interactions[d11_inter.get_category()][d11_inter.enrichment_status_tag_pair].append(d11_inter)
        self.interactions['ALL'][d11_inter.enrichment_status_tag_pair].append(d11_inter)

    def get_all_pairwise_differences_of_interaction_distances(self, i_cat, e_cat):

        d_inter_list = self.interactions[i_cat][e_cat]
        d_inter_list_len = len(d_inter_list)

        # Get list of all pairwise interaction distances
        pairwise_i_dist_diffs = []
        for i in range(0, d_inter_list_len):
            dist_a = d_inter_list[i].i_dist
            for j in range(i + 1, d_inter_list_len):
                dist_b = d_inter_list[j].i_dist
                diff = abs(dist_a - dist_b)
                pairwise_i_dist_diffs.append(diff)

        return pairwise_i_dist_diffs

    def get_interaction_number(self, i_cat, e_cat):
        return len(self.interactions[i_cat][e_cat])

    def get_read_pair_number(self, i_cat, e_cat):

        # Get list of interactions
        d_inter_list = self.interactions[i_cat][e_cat]

        # Sum up read pairs numbers of interactions
        rp_total_sum = 0
        for d_inter in d_inter_list:
            rp_total_sum += d_inter.rp_total

        return rp_total_sum

    def get_median_interaction_distance(self, i_cat, e_cat):

        # Get list of interactions
        d_inter_list = self.interactions[i_cat][e_cat]

        # Collect interaction distances
        distances = []
        for d_inter in d_inter_list:
            distances.append(d_inter.i_dist)

        if 0 < len(distances):
            mean_interaction_distance = int(np.median(distances))
        else:
            mean_interaction_distance = 0

        return mean_interaction_distance

    def n_total_interactions(self):
        return len(self.interactions['ALL']['NE']) + len(self.interactions['ALL']['EN'])

    def n_di_interactions(self):
        return len(self.interactions['DI']['NE']) + len(self.interactions['DI']['EN'])

    def n_uir_interactions(self):
        return len(self.interactions['UIR']['NE']) + len(self.interactions['UIR']['EN'])

    def n_di_ne_interactions(self):
        return len(self.interactions['DI']['NE'])

    def n_uir_ne_interactions(self):
        return len(self.interactions['UIR']['NE'])

    def get_curb_number(self, i_cat, e_cat, curb_size: int = 270500, curb_size_error_margin: int = 10000):

        if i_cat not in self.curb_nums:
            self.curb_nums[i_cat] = dict()

        if e_cat not in self.curb_nums[i_cat]:
            self.curb_nums[i_cat][e_cat] = dict()

        # Get interactions
        d_inter_list = self.interactions[i_cat][e_cat]
        d_inter_list_len = len(d_inter_list)

        # Init counter
        curb_num = 0

        # Calculate all pairwise differences of interaction distances
        for i in range(0, d_inter_list_len):
            dist_a = d_inter_list[i].i_dist
            for j in range(i + 1, d_inter_list_len):
                dist_b = d_inter_list[j].i_dist
                diff = abs(dist_a - dist_b)
                # Count differences that are a multiple of the curb size
                remainder = diff % curb_size
                if curb_size - curb_size_error_margin < diff and\
                        (remainder < curb_size_error_margin or curb_size - remainder < curb_size_error_margin):
                    curb_num += 1
                    #print(str(remainder) + '\t' + str(diff) + '\t' + str(curb_num))
        self.curb_nums[i_cat][e_cat] = curb_num
        return curb_num
