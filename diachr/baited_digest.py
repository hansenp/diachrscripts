from .diachromatic_interaction import DiachromaticInteraction11
import numpy as np


class BaitedDigest:
    """
    Objects of this class are used to group interactions interactions that end in the same bait. Interactions are broken
    down by category and direction as seen from the bait. For the interaction categories unbalanced (DI), balanced
    reference (UIR), balanced (UI) and all interactions (ALL) there are two lists of interactions each, one with
    interactions from the bait to the left (NE) and one with interactions from the bait to the right (EN).
    """

    # Dictionary with lists of interactions at this baited digests
    interactions = None

    curb_nums = None

    def __init__(self):

        self.interactions = {
            'DIX': {
                'NE': [],
                'EN': []},
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

    def get_rp_num_list(self, i_cat, e_cat):

        # Get list of interactions
        d_inter_list = self.interactions[i_cat][e_cat]
        rp_num_list = []
        for d_inter in d_inter_list:
            rp_num_list.append(d_inter.rp_total)
        return rp_num_list

    def get_i_dist_list(self, i_cat, e_cat):

        # Get list of interactions
        d_inter_list = self.interactions[i_cat][e_cat]
        i_dist_list = []
        for d_inter in d_inter_list:
            i_dist_list.append(d_inter.i_dist)
        return i_dist_list

    def get_dig_len_list(self, i_cat, return_shorter_len: bool = True):
        """
        Returns a list of digest lengths of all interactions in a category, either always the shorter or longer length.

        :param i_cat: Interaction category.
        :param return_shorter_len: If true, return the length  of the shorter digest for each interaction, otherwise the
        length of the longer digest.
        :return: List of digest lengths of all interactions.
        """

        # Get list of interactions
        d_inter_list = self.interactions[i_cat]['NE'] + self.interactions[i_cat]['EN']

        # Init list that will be returned
        dig_len_list = []

        # Iterate interactions
        for d_inter in d_inter_list:

            # Get digests lengths
            d1_len = d_inter._toA - d_inter._fromA
            d2_len = d_inter._toB - d_inter._fromB

            # Add length of shorter or longer digest to list
            if return_shorter_len:
                if d1_len < d2_len:
                    dig_len_list.append(d1_len)
                else:
                    dig_len_list.append(d2_len)
            else:
                if d1_len < d2_len:
                    dig_len_list.append(d2_len)
                else:
                    dig_len_list.append(d1_len)

        return dig_len_list

    def get_median_read_pair_number(self, i_cat, e_cat):

        # Get list of interactions
        d_inter_list = self.interactions[i_cat][e_cat]

        # Collect read pair numbers of interactions
        rp_list = []
        for d_inter in d_inter_list:
            rp_list.append(d_inter.rp_total)

        # Get median read pair number
        if 0 < len(rp_list):
            return np.median(rp_list)
        else:
            return 0.0

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

    def get_interactions_sorted_by_dist(self, i_cat, e_cat):

        distances = []
        for d11_inter in self.interactions[i_cat][e_cat]:
            distances.append(d11_inter.i_dist)
        sorted_distance_idx = np.argsort(distances)
        sorted_interaction_list = []
        for idx in sorted_distance_idx:
            sorted_interaction_list.append(self.interactions[i_cat][e_cat][idx])
        return sorted_interaction_list
