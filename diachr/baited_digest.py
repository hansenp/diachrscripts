from .diachromatic_interaction import DiachromaticInteraction11

class BaitedDigest:
    """
    Objects of this class are used to group interactions interactions that end in the same bait. Interactions are broken
    down by category and direction as seen from the bait. For the interaction categories directed (DI), undirected
    reference (UIR), undirected (UI) and all interactions (ALL) there are two lists of interactions each, one with
    interactions from the bait to the left (NE) and one with interactions from the bait to the right (EN).
    """

    # Dictionary with lists of interactions ending in baited digest
    interactions = None

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

    def add_interaction(self, d11_inter: DiachromaticInteraction11):
        self.interactions[d11_inter.get_category()][d11_inter.enrichment_status_tag_pair].append(d11_inter)
        self.interactions['ALL'][d11_inter.enrichment_status_tag_pair].append(d11_inter)

    def get_all_pairwise_differences_of_interaction_distances(self, i_cat, dir):

        d_inter_list = self.interactions[i_cat][dir]
        d_inter_list_len = len(d_inter_list)

        # Get list of all pairwise interaction distances
        pairwise_i_dist_diffs = []
        for i in range(0, d_inter_list_len):
            dist_a = d_inter_list[i].i_dist
            for j in range(i + 1, d_inter_list_len):
                dist_b = d_inter_list[j].i_dist
                pairwise_i_dist_diffs.append(abs(dist_a - dist_b))

        return pairwise_i_dist_diffs

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
