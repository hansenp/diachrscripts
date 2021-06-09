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

    def add_di_simple_twisted_categories(self):
        """
        This function divides the directed interactions into the categories simple and twisted and creates two new
        corresponding dictionaries.
        """

        self.interactions['DI_S'] = {'NE': [], 'EN': []}
        self.interactions['DI_T'] = {'NE': [], 'EN': []}

        # Get lists of interactions
        d_inter_list = self.interactions['DI']['NE'] + self.interactions['DI']['EN']

        for d_inter in d_inter_list:
            if d_inter.n_twisted < d_inter.n_simple:
                self.interactions['DI_S'][d_inter.enrichment_status_tag_pair].append(d_inter)
            else:
                self.interactions['DI_T'][d_inter.enrichment_status_tag_pair].append(d_inter)

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

    def select_undirected_reference_interactions(self, lower_q: float = 0.25, upper_q: float = 0.75):
        """
        This function selects undirected reference interactions at this bait that, as compared to the directed at this
        bait, have a similar distribution of read pairs.

        :param lower_q: Lower quantile
        :param upper_q: Upper quantile
        """

        # Create list of read pair numbers from directed interactions
        di_rp_list = []
        for d11_inter in self.interactions['DI']['NE'] + self.interactions['DI']['EN']:
            di_rp_list.append(d11_inter.rp_total)

        # Determine lower and upper quantile for read pair numbers from directed interactions
        if 0 < len(di_rp_list):
            q1 = np.quantile(di_rp_list, lower_q)
            q3 = np.quantile(di_rp_list, upper_q)
        else:
            # If there are no directed interactions, we do not select reference interactions
            q1 = 0.0
            q3 = 0.0

        # Concatenate all lists with undirected interactions at this bait (without assumptions about NE-EN distribution)
        ui_list = self.interactions['UI']['NE'] + self.interactions['UI']['EN'] + \
                  self.interactions['UIR']['NE'] + self.interactions['UIR']['EN']

        # Select reference interactions with read pair numbers between lower and upper quantile
        uir_rp_dict = {
            'NE': [],
            'EN': []
        }
        uir_dict = {
            'NE': [],
            'EN': []
        }
        for d11_inter in ui_list:
            if q1 <= d11_inter.rp_total <= q3:
                uir_dict[d11_inter.enrichment_status_tag_pair].append(d11_inter)
                uir_rp_dict[d11_inter.enrichment_status_tag_pair].append(d11_inter.rp_total)

        # For now, assign the reference interactions to the variables for all interactions so that we can compare the
        # old and new reference
        self.interactions['ALL']['NE'] = uir_dict['NE']
        self.interactions['ALL']['EN'] = uir_dict['EN']

        print("-------")
        print("Number of DI: " + str(len(di_rp_list)))
        print("DI read pair numbers: " + str(sorted(di_rp_list)))
        print("Upper quantile: " + str(q1))
        print("Lower quantile: " + str(q3))
        if 0 < len(di_rp_list):
            print("Median read pair number DI: " + str(np.median(di_rp_list)))
        else:
            print('NA')
        if 0 < len(uir_rp_dict['NE']):
            print("Median read pair number UIR-NE: " + str(np.median(uir_rp_dict['NE'])))
        else:
            print('NA')
        if 0 < len(uir_rp_dict['EN']):
            print("Median read pair number UIR-EN: " + str(np.median(uir_rp_dict['EN'])))
        else:
            print('NA')

    def get_interactions_sorted_by_dist(self, i_cat, e_cat):

        distances = []
        for d11_inter in self.interactions[i_cat][e_cat]:
            distances.append(d11_inter.i_dist)
        sorted_distance_idx = np.argsort(distances)
        sorted_interaction_list = []
        for idx in sorted_distance_idx:
            sorted_interaction_list.append(self.interactions[i_cat][e_cat][idx])
        return sorted_interaction_list
