from .baited_digest import BaitedDigest
from .diachromatic_interaction import DiachromaticInteraction11
from .diachromatic_interaction_set import DiachromaticInteractionSet
from collections import defaultdict
import matplotlib.pyplot as plt

class BaitedDigestSet:

    def __init__(self):

        # Dictionary that contains all BaitedDigest objects
        self._baited_digest_dict = defaultdict()

        # Dictionary with infos about ingestion of interactions
        self._ingest_interaction_set_info_dict = {
            'TOTAL_INTERACTIONS_READ': 0,
            'DI': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'UIR': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'UI': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'ALL': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'BAITED_DIGESTS': 0
        }

    def ingest_interaction_set(self, d11_inter_set: DiachromaticInteractionSet, verbose: bool = False):

        if verbose:
            print("[INFO] Reading interactions and group them according to chromosomes and baited digests ...")

        for d11_inter in d11_inter_set.interaction_list:
            self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] += 1
            if verbose:
                if self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] % 1000000 == 0:
                    print("\t[INFO] Read " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + " interactions ...")

            # Count interaction type
            self._ingest_interaction_set_info_dict[d11_inter.get_category()][d11_inter.enrichment_status_tag_pair] += 1
            self._ingest_interaction_set_info_dict['ALL'][d11_inter.enrichment_status_tag_pair] += 1

            # Skip NN and EE interactions
            if d11_inter.enrichment_status_tag_pair == 'NN' or d11_inter.enrichment_status_tag_pair == 'EE':
                continue

            # Get key from interaction
            if d11_inter.enrichment_status_tag_pair == 'EN':
                enriched_digest_coords = d11_inter.chrA + '\t' + str(d11_inter.fromA) + '\t' + str(d11_inter.toA)
            if d11_inter.enrichment_status_tag_pair == 'NE':
                enriched_digest_coords = d11_inter.chrB + '\t' + str(d11_inter.fromB) + '\t' + str(d11_inter.toB)

            # Create a new dictionary if this is the first interaction seen on this chromosome
            if d11_inter.chrA not in self._baited_digest_dict:
                self._baited_digest_dict[d11_inter.chrA] = defaultdict()

            # Create new BaitedDigest object if it doesn't already exists
            if enriched_digest_coords not in self._baited_digest_dict[d11_inter.chrA].keys():
                self._baited_digest_dict[d11_inter.chrA][enriched_digest_coords] = BaitedDigest()
                self._ingest_interaction_set_info_dict['BAITED_DIGESTS'] += 1

            # Add interaction BaitedDigest object
            self._baited_digest_dict[d11_inter.chrA][enriched_digest_coords].add_interaction(d11_inter)

        if verbose:
            print("\t[INFO] Total number of interactions read: " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']))
            print("\t[INFO] Total number of baited digests: " + "{:,}".format(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']))
            print("[INFO] ... done.")

    def get_ingest_interaction_set_info_dict(self):
        """
        :return: Dictionary that contains information about parsed interaction data.
        """
        return self._ingest_interaction_set_info_dict

    def get_ingest_interaction_set_info_report(self):

        report = "[INFO] Report on ingestion of interactions:" + '\n'
        report += "\t[INFO] Total number of interactions read: " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + '\n'
        discarded = self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] - \
                    self._ingest_interaction_set_info_dict['ALL']['NE'] - \
                    self._ingest_interaction_set_info_dict['ALL']['EN']
        report += "\t[INFO] Discarded NN and EE interactions: " + "{:,}".format(discarded) + '\n'
        ingested = self._ingest_interaction_set_info_dict['ALL']['NE'] + \
                   self._ingest_interaction_set_info_dict['ALL']['EN']
        report += "\t[INFO] Total number of ingested NE and EN interactions: " + "{:,}".format(ingested) + '\n'
        report += "\t[INFO] Broken down by interaction category and enrichment status: " + '\n'
        report += "\t\t[INFO] DI: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['DI']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['DI']['EN']) + '\n'
        report += "\t\t[INFO] UIR: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['UIR']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['UIR']['EN']) + '\n'
        report += "\t\t[INFO] UI: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['UI']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['UI']['EN']) + '\n'
        report += "\t\t[INFO] ALL: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['ALL']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['ALL']['EN']) + '\n'
        report += "\t[INFO] Total number of baited digests: " + "{:,}".format(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_ingest_interaction_set_table_row(self):

        # Header row
        table_row = ":TR_INGESTION:" + '\t'

        table_row += "TOTAL_INTERACTIONS_READ" + '\t'
        table_row += "DISCARDED" + '\t'
        table_row += "INGESTED" + '\t'

        table_row += "DI_NE" + '\t'
        table_row += "DI_EN" + '\t'

        table_row += "UIR_NE" + '\t'
        table_row += "UIR_EN" + '\t'

        table_row += "UI_NE" + '\t'
        table_row += "UI_EN" + '\t'

        table_row += "ALL_NE" + '\t'
        table_row += "ALL_EN" + '\t'

        table_row += "BAITED_DIGESTS" + '\n'

        # Row with values
        table_row += ":TR_INGESTION:" + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + '\t'
        discarded = self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] - \
                    self._ingest_interaction_set_info_dict['ALL']['NE'] - \
                    self._ingest_interaction_set_info_dict['ALL']['EN']
        table_row += str(discarded) + '\t'
        ingested = self._ingest_interaction_set_info_dict['ALL']['NE'] + \
                   self._ingest_interaction_set_info_dict['ALL']['EN']
        table_row += str(ingested) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['DI']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['DI']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['UIR']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['UIR']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['UI']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['UI']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['ALL']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['ALL']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']) + '\n'

        return table_row

    def get_baited_digest_keys_sorted_by_sta_pos(self):
        """
        :return: A dictionary with lists of baited digests keys sorted by starting coordinates, one for each chromosome.
        """

        # Create a dictionary with starting coordinates as key and digest coordinates as values
        sta_keys_dict = defaultdict()
        for chr in self._baited_digest_dict.keys():
            sta_keys_dict[chr] = defaultdict()
            for enriched_digest_coords in self._baited_digest_dict[chr].keys():
                chr_bd, sta_bd, end_bd = enriched_digest_coords.split('\t')
                sta_keys_dict[chr][sta_bd] = enriched_digest_coords

        # Create sorted lists of baited digests keys for each chromosome
        sorted_digest_key_lists = defaultdict()
        for chr in self._baited_digest_dict.keys():
            sorted_digest_key_list = []
            for sta_key in sorted(sta_keys_dict[chr].keys()):
                sorted_digest_key_list.append(sta_keys_dict[chr][sta_key])
            sorted_digest_key_lists[chr] = sorted_digest_key_list

        return sorted_digest_key_lists

    def proportion_of_directed_interactions_on_individual_chromosomes(self):

        for chr in self._baited_digest_dict.keys():
            print(len(self._baited_digest_dict[chr]))
            n_total_interactions = 0
            n_di_interactions = 0
            n_uir_interactions = 0
            n_di_ne_interactions = 0
            n_uir_ne_interactions = 0
            for key, baited_digest in self._baited_digest_dict[chr].items():
                n_total_interactions += baited_digest.n_total_interactions()
                n_di_interactions += baited_digest.n_di_interactions()
                n_uir_interactions += baited_digest.n_uir_interactions()
                n_di_ne_interactions += baited_digest.n_di_ne_interactions()
                n_uir_ne_interactions += baited_digest.n_uir_ne_interactions()
            proportion_di = round(n_di_interactions/float(n_total_interactions), 2)
            proportion_uir = round(n_uir_interactions / float(n_total_interactions), 2)
            di_ne_proportion = round(n_di_ne_interactions/float(n_di_interactions), 2)
            uir_ne_proportion = round(n_uir_ne_interactions/float(n_uir_interactions), 2)
            print(chr + '\t' +
                  str(n_total_interactions) + '\t' +
                  str(n_di_interactions) + '\t' +
                  str(proportion_di) + '\t' +
                  str(proportion_uir) + '\t' +
                  str(di_ne_proportion) + '\t' +
                  str(uir_ne_proportion)
                  )

    def get_pairwise_interaction_distances_at_baits(self, chromosomes: [str] = None):

        # Prepare data structure for results
        i_cats = ['DI', 'UIR', 'UI', 'ALL']
        e_cats = ['NE', 'EN']
        pid = dict()
        for i_cat in i_cats:
            pid[i_cat] = dict()
            for e_cat in e_cats:
                pid[i_cat][e_cat] = []
        pid['CHROMOSOMES'] = []

        # Combine lists of pairwise distances from all baits
        for chr in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chr not in chromosomes:
                    continue

            pid['CHROMOSOMES'].append(chr)
            for key, baited_digest in self._baited_digest_dict[chr].items():
                for i_cat in i_cats:
                    for e_cat in e_cats:
                        pid[i_cat][e_cat].extend(baited_digest.get_all_pairwise_differences_of_interaction_distances(i_cat, e_cat))

        return pid

    def get_pairwise_interaction_distances_at_baits_histograms(self,
                                                               pid_dict: dict = None,
                                                               description: str = "DESCRIPTION",
                                                               pdf_file_name: str = "pid_histograms.pdf"):

        i_cats = ['DI', 'UIR', 'UI', 'ALL']

        fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(11, 19),gridspec_kw={'height_ratios': [0.5, 1, 1, 1, 1, 1, 1, 1]})

        ax[0][0].plot()
        ax[0][0].spines['left'].set_color('white')
        ax[0][0].spines['right'].set_color('white')
        ax[0][0].spines['top'].set_color('white')
        ax[0][0].spines['bottom'].set_color('white')
        ax[0][0].tick_params(axis='x', colors='white')
        ax[0][0].tick_params(axis='y', colors='white')
        ax[0][1].plot()
        ax[0][1].spines['left'].set_color('white')
        ax[0][1].spines['right'].set_color('white')
        ax[0][1].spines['top'].set_color('white')
        ax[0][1].spines['bottom'].set_color('white')
        ax[0][1].tick_params(axis='x', colors='white')
        ax[0][1].tick_params(axis='y', colors='white')

        fig.text(0.015, 0.9825, 'Pairwise differences of interaction distances at baits', fontsize=18, fontweight='bold')
        fig.text(0.030, 0.9660, 'Description: ' + description, fontsize=12)
        fig.text(0.030, 0.9525, 'For chromosomes:', fontsize=12)
        fig.text(0.045, 0.9425, str(pid_dict['CHROMOSOMES']), fontsize=8)

        # Create histograms
        x_ticks = [0, 500000, 1000000, 1500000, 2000000]
        bin_width = 20000
        i_cat_colors = {'DI': (255 / 255, 163 / 255, 0 / 255, 1), 'UIR': (171 / 255, 215 / 255, 230 / 255, 1),
                        'UI': (210 / 255, 210 / 255, 210 / 255, 1), 'ALL': 'black'}
        i_cat_names = {'DI': 'Directed', 'UIR': 'Undirected reference', 'UI': 'Undirected', 'ALL': 'All'}

        # Create histograms for DI, UIR, UI and ALL
        for i in [0, 1, 2, 3]:

            # Prepare bins
            x_max = max(pid_dict[i_cats[i]]['NE'] + pid_dict[i_cats[i]]['EN'])
            bins = range(0, x_max + bin_width, bin_width)

            # Create histogram for NE
            counts_1, bins, patches = ax[i+1][0].hist(
                pid_dict[i_cats[i]]['NE'],
                bins=bins, density=False,
                facecolor=i_cat_colors[i_cats[i]],
                edgecolor="black", alpha=0.75)
            ax[i+1][0].set_title(i_cat_names[i_cats[i]] + ' to the left (n=' + "{:,}".format(len(pid_dict[i_cats[i]]['NE'])) + ')', loc='left')
            ax[i+1][0].set_xlabel('Pairwise differences of interaction distances')
            ax[i+1][0].set_ylabel('Frequency')
            ax[i+1][0].set_xticks(x_ticks)
            ax[i+1][0].set_xlim(0, 2000000)

            # Create histogram for NE
            counts_2, bins, patches = ax[i+1][1].hist(
                pid_dict[i_cats[i]]['EN'],
                bins=bins, density=False,
                facecolor=i_cat_colors[i_cats[i]],
                edgecolor="black", alpha=0.75)
            ax[i+1][1].set_title(i_cat_names[i_cats[i]] + ' to the right (n=' + "{:,}".format((len(pid_dict[i_cats[i]]['EN']))) + ')', loc='left')
            ax[i+1][1].set_xlabel('Pairwise differences of interaction distances')
            ax[i+1][1].set_ylabel('Frequency')
            ax[i+1][1].set_xticks(x_ticks)
            ax[i+1][1].set_xlim(0, 2000000)

            # Make y-axes comparable
            ymax = max(max(counts_1), max(counts_2))
            ax[i+1][0].set_ylim(0, ymax)
            ax[i+1][1].set_ylim(0, ymax)

        # Draw vertical lines for UIR - NE
        ax[1+1][0].axvline(1 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(2 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(3 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(4 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(5 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(6 * 270500, color='blue', linewidth=0.5, zorder=0)
        ax[1+1][0].axvline(7 * 270500, color='blue', linewidth=0.5, zorder=0)

        first_count_bin = 5
        for i in [5, 6, 7]:

            # Prepare bins
            x_max = max(pid_dict['UIR']['NE'] + pid_dict['UIR']['EN'])
            bins = range(0, x_max + bin_width, bin_width)

            # Create histogram for NE
            counts_1, bins, patches = ax[i][0].hist(
                pid_dict['UIR']['NE'],
                bins=bins, density=False,
                facecolor=i_cat_colors['UIR'],
                edgecolor="black", alpha=0.75)
            ax[i][0].set_title(i_cat_names['UIR'] + ' to the left (n=' + "{:,}".format(len(pid_dict['UIR']['NE'])) + ')',
                                   loc='left')
            ax[i][0].set_xlabel('Pairwise differences of interaction distances')
            ax[i][0].set_ylabel('Frequency')
            ax[i][0].set_xticks(x_ticks)
            ax[i][0].set_xlim(0, 2000000)

            # Create histogram for NE
            counts_2, bins, patches = ax[i][1].hist(
                pid_dict['UIR']['EN'],
                bins=bins, density=False,
                facecolor=i_cat_colors['UIR'],
                edgecolor="black", alpha=0.75)
            ax[i][1].set_title(i_cat_names['UIR'] + ' to the right (n=' + "{:,}".format(len(pid_dict['UIR']['EN'])) + ')',
                                   loc='left')
            ax[i][1].set_xlabel('Pairwise differences of interaction distances')
            ax[i][1].set_ylabel('Frequency')
            ax[i][1].set_xticks(x_ticks)
            ax[i][1].set_xlim(0, 2000000)

            # Make y-axes comparable
            ymax = max(counts_1[first_count_bin:])
            ax[i][0].set_ylim(0, ymax)
            ax[i][1].set_ylim(0, ymax)
            first_count_bin += 14

            # Draw vertical lines for UIR - NE
            ax[i][0].axvline(1 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(2 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(3 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(4 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(5 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(6 * 270500, color='blue', linewidth=0.5, zorder=0)
            ax[i][0].axvline(7 * 270500, color='blue', linewidth=0.5, zorder=0)

        # Save figure
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig




