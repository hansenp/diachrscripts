from .baited_digest import BaitedDigest
from .diachromatic_interaction_set import DiachromaticInteractionSet
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np


class BaitedDigestSet:

    def __init__(self):

        # Dictionary that contains ingested interaction grouped by chromosomes and baits
        self._baited_digest_dict = defaultdict()

        # Dictionary with information about ingestion of interactions
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
        """
        Ingests interactions from a DiachromaticInteractionSet and groups them by chromosomes and baits.
        :param d11_inter_set: DiachromaticInteractionSet with DiachromaticInteraction11 interactions
        :param verbose: If true, progress information will be displayed
        :return: Dictionary containing information about ingested  interactions
        """

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

        return self._ingest_interaction_set_info_dict

    def get_ingest_interaction_set_info_report(self):
        """
        :return: Formatted string with information about ingestion of interactions for output on the screen or in report
         files.
        """

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
        """
        :return: Table row with information about ingestion of interactions
        """

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

    def get_pairwise_differences_of_interaction_distances_at_baits(self, chromosomes: [str] = None, verbose: bool = False):
        """
        This function performs the analysis with all pairwise differences of interaction distances at baits. The
        pairwise distances are determined within the individual BaitedDigest objects and returned as lists. In this,
        function, the pairwise distances from all baits are combined.
        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21']
        :return: A dictionary containing the results and a list of chromosomes that were taken into account.
        """

        if verbose:
            print("[INFO] Getting all pairwise differences of interaction distances at baits ...")

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

            if verbose:
                print("\t[INFO] Processing chromosome " + chr + " ...")

            pid['CHROMOSOMES'].append(chr)
            for key, baited_digest in self._baited_digest_dict[chr].items():
                for i_cat in i_cats:
                    for e_cat in e_cats:
                        pid[i_cat][e_cat].extend(baited_digest.get_all_pairwise_differences_of_interaction_distances(i_cat, e_cat))

        if verbose:
            print("[INFO] ... done.")

        return pid

    def get_pairwise_differences_of_interaction_distances_at_baits_histograms(self,
                                                               pid_dict: dict = None,
                                                               description: str = "DESCRIPTION",
                                                               pdf_file_name: str = "pid_histograms.pdf"):
        """
        This function creates the histograms for the pairwise interaction distances at baits for all interaction
        categories and enrichment states.
        :param pid_dict: Dictionary that was created with the function 'get_pairwise_interaction_distances_at_baits()'
        :param description: Brief description that is shown in the plot above the historgrams
        :param pdf_file_name: Name of the PDF file that will be created.
        :return: A matplotlib 'Figure' object that can be displayed in Jupyter notebooks
        """

        # Interaction categories
        i_cats = ['DI', 'UIR', 'UI', 'ALL']

        # Prepare grid for individual plots
        fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(11, 19),gridspec_kw={'height_ratios': [0.5, 1, 1, 1, 1, 1, 1, 1]})

        # Add header section with description and chromosomes that were taken into account
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

        # Set variables that all historgrams have in common
        x_lim = 2000000
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
            ax[i+1][0].set_xlim(0, x_lim)

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
            ax[i+1][1].set_xlim(0, x_lim)

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

        # Add additional histograms for UIR with smaller limits on the y-axis
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
            ax[i][0].set_xlim(0, x_lim)

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
            ax[i][1].set_xlim(0, x_lim)

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

        # Save and return figure
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    def get_interaction_numbers_at_baits(self, chromosomes: [str] = None):
        """
        This function performs the analysis with interaction at baits. The
        interaction numbers are determined within the individual BaitedDigest objects as numbers. In this,
        function, the interaction numbers from all baits are combined.
        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21']
        :return: A dictionary containing the results and a list of chromosomes that were taken into account.
        """

        # Prepare data structure for results
        i_cats = ['DI', 'UIR', 'UI', 'ALL']
        e_cats = ['NE', 'EN']
        i_nums = dict()
        for i_cat in i_cats:
            i_nums[i_cat] = dict()
            for e_cat in e_cats:
                i_nums[i_cat][e_cat] = []
        i_nums['CHROMOSOMES'] = []

        # Combine lists of pairwise distances from all baits
        for chr in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chr not in chromosomes:
                    continue

            i_nums['CHROMOSOMES'].append(chr)
            for key, baited_digest in self._baited_digest_dict[chr].items():
                for i_cat in i_cats:
                    for e_cat in e_cats:
                        i_nums[i_cat][e_cat].append(baited_digest.get_interaction_number(i_cat, e_cat))

        return i_nums

    def get_interaction_numbers_at_baits_histograms(self,
                                                               i_num_dict: dict = None,
                                                               description: str = "DESCRIPTION",
                                                               pdf_file_name: str = "i_num_histograms.pdf"):
        """
        This function creates the histograms for interaction numbers at baits for all interaction
        categories and enrichment states.
        :param pid_dict: Dictionary that was created with the function 'get_interaction_numbers_at_baits()'
        :param description: Brief description that is shown in the plot above the historgrams
        :param pdf_file_name: Name of the PDF file that will be created.
        :return: A matplotlib 'Figure' object that can be displayed in Jupyter notebooks
        """

        # Interaction categories
        i_cats = ['DI', 'UIR', 'UI', 'ALL']

        # Prepare grid for individual plots
        fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(11, 19),gridspec_kw={'height_ratios': [0.5, 1, 1, 1, 1, 1, 1, 1]})

        # Add header section with description and chromosomes that were taken into account
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
        fig.text(0.015, 0.9825, 'Interaction numbers at baits', fontsize=18, fontweight='bold')
        fig.text(0.030, 0.9660, 'Description: ' + description, fontsize=12)
        fig.text(0.030, 0.9525, 'For chromosomes:', fontsize=12)
        fig.text(0.045, 0.9425, str(i_num_dict['CHROMOSOMES']), fontsize=8)

        # Set variables that all historgrams have in common
        x_lim = 200
        x_ticks = [0, 50, 100, 150, 200]
        bin_width = 1
        i_cat_colors = {'DI': (255 / 255, 163 / 255, 0 / 255, 1), 'UIR': (171 / 255, 215 / 255, 230 / 255, 1),
                        'UI': (210 / 255, 210 / 255, 210 / 255, 1), 'ALL': 'black'}
        i_cat_names = {'DI': 'Directed', 'UIR': 'Undirected reference', 'UI': 'Undirected', 'ALL': 'All'}

        # Create histograms for DI, UIR, UI and ALL
        for i in [0, 1, 2, 3]:

            # Prepare bins
            x_max = max(i_num_dict[i_cats[i]]['NE'] + i_num_dict[i_cats[i]]['EN'])
            bins = range(0, x_max + bin_width, bin_width)

            # Create histogram for NE
            counts_1, bins, patches = ax[i+1][0].hist(
                i_num_dict[i_cats[i]]['NE'],
                bins=bins, density=False,
                facecolor=i_cat_colors[i_cats[i]],
                edgecolor="black", alpha=0.75)
            ax[i+1][0].set_title(i_cat_names[i_cats[i]] + ' to the left (n=' + "{:,}".format(sum(i_num_dict[i_cats[i]]['NE'])) + ')', loc='left')
            ax[i+1][0].set_xlabel('Interaction numbers')
            ax[i+1][0].set_ylabel('Frequency')
            ax[i+1][0].set_xticks(x_ticks)
            ax[i+1][0].set_xlim(0, x_lim)

            # Create histogram for NE
            counts_2, bins, patches = ax[i+1][1].hist(
                i_num_dict[i_cats[i]]['EN'],
                bins=bins, density=False,
                facecolor=i_cat_colors[i_cats[i]],
                edgecolor="black", alpha=0.75)
            ax[i+1][1].set_title(i_cat_names[i_cats[i]] + ' to the right (n=' + "{:,}".format((sum(i_num_dict[i_cats[i]]['EN']))) + ')', loc='left')
            ax[i+1][1].set_xlabel('Interaction numbers')
            ax[i+1][1].set_ylabel('Frequency')
            ax[i+1][1].set_xticks(x_ticks)
            ax[i+1][1].set_xlim(0, x_lim)

            # Make y-axes comparable
            ymax = max(max(counts_1), max(counts_2))
            ax[i+1][0].set_ylim(0, ymax)
            ax[i+1][1].set_ylim(0, ymax)

        # Add additional histograms for UIR with smaller limits on the y-axis
        divisor = 1
        for i in [5, 6, 7]:

            # Prepare bins
            x_max = max(i_num_dict['UIR']['NE'] + i_num_dict['UIR']['EN'])
            bins = range(0, x_max + bin_width, bin_width)

            # Create histogram for NE
            counts_1, bins, patches = ax[i][0].hist(
                i_num_dict['UIR']['NE'],
                bins=bins, density=False,
                facecolor=i_cat_colors['UIR'],
                edgecolor="black", alpha=0.75)
            ax[i][0].set_title(i_cat_names['UIR'] + ' to the left (n=' + "{:,}".format(sum(i_num_dict['UIR']['NE'])) + ')',
                                   loc='left')
            ax[i][0].set_xlabel('Pairwise differences of interaction distances')
            ax[i][0].set_ylabel('Frequency')
            ax[i][0].set_xticks(x_ticks)
            ax[i][0].set_xlim(0, x_lim)

            # Create histogram for NE
            counts_2, bins, patches = ax[i][1].hist(
                i_num_dict['UIR']['EN'],
                bins=bins, density=False,
                facecolor=i_cat_colors['UIR'],
                edgecolor="black", alpha=0.75)
            ax[i][1].set_title(i_cat_names['UIR'] + ' to the right (n=' + "{:,}".format(sum(i_num_dict['UIR']['EN'])) + ')',
                                   loc='left')
            ax[i][1].set_xlabel('Pairwise differences of interaction distances')
            ax[i][1].set_ylabel('Frequency')
            ax[i][1].set_xticks(x_ticks)
            ax[i][1].set_xlim(0, x_lim)

            # Make y-axes comparable
            ymax = max(counts_1)/divisor
            divisor *= 8
            ax[i][0].set_ylim(0, ymax)
            ax[i][1].set_ylim(0, ymax)

        # Save and return figure
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    def get_interaction_number_pairs_at_baits(self, chromosomes: [str] = None):

        # Prepare data structure for results
        i_cats = ['DI', 'UIR', 'UI', 'ALL']
        e_cats = ['NE', 'EN']
        i_num_pairs = dict()
        for i_cat in i_cats:
            i_num_pairs[i_cat] = dict()
            for e_cat in e_cats:
                i_num_pairs[i_cat][e_cat] = []
        i_num_pairs['CHROMOSOMES'] = []

        # Combine lists of pairwise distances from all baits
        for chr in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chr not in chromosomes:
                    continue

            i_num_pairs['CHROMOSOMES'].append(chr)
            for key, baited_digest in self._baited_digest_dict[chr].items():
                for i_cat in i_cats:
                    i_num_pairs[i_cat]['NE'].append(baited_digest.get_interaction_number(i_cat, 'NE'))
                    i_num_pairs[i_cat]['EN'].append(baited_digest.get_interaction_number(i_cat, 'EN'))

        return i_num_pairs

    def get_read_pair_number_pairs_at_baits(self, chromosomes: [str] = None):

        # Prepare data structure for results
        i_cats = ['DI', 'UIR', 'UI', 'ALL']
        e_cats = ['NE', 'EN']
        rp_num_pairs = dict()
        for i_cat in i_cats:
            rp_num_pairs[i_cat] = dict()
            for e_cat in e_cats:
                rp_num_pairs[i_cat][e_cat] = []
        rp_num_pairs['CHROMOSOMES'] = []

        # Combine lists of pairwise distances from all baits
        for chr in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chr not in chromosomes:
                    continue

            rp_num_pairs['CHROMOSOMES'].append(chr)
            for key, baited_digest in self._baited_digest_dict[chr].items():
                for i_cat in i_cats:
                    rp_num_pairs[i_cat]['NE'].append(baited_digest.get_read_pair_number(i_cat, 'NE'))
                    rp_num_pairs[i_cat]['EN'].append(baited_digest.get_read_pair_number(i_cat, 'EN'))

        return rp_num_pairs




#############

    def make_ticks(self, max_val: int=100):
        if max_val < 100:
            tick_dist = 20
        elif max_val < 500:
            tick_dist = 50
        elif max_val < 1000:
            tick_dist = 200
        elif max_val < 1500:
            tick_dist = 500
        elif max_val < 5000:
            tick_dist = 1000
        elif max_val < 10000:
            tick_dist = 2000
        elif max_val < 25000:
            tick_dist = 5000
        else:
            tick_dist = 10000
        ticks = range(0, max_val + tick_dist, tick_dist)
        if max_val < 1000:
            tick_lables = ticks
        else:
            tick_lables = []
            for tick in ticks:
                tick_lables.append(str(tick/1000) + 'k')

        return ticks, tick_lables

    def custom_pair_plot(self,
                         i_cat='I_CAT',
                         en=None,
                         ne=None,
                         color='black',
                         bin_size=1,
                         xy_max=None,
                         ax_hx=None,
                         ax_hy=None,
                         ax_s=None,
                         **plt_kwargs):
        if ax_hx is None:
            ax_hx = plt.gca()
        if ax_hy is None:
            ax_hy = plt.gca()
        if ax_s is None:
            ax_s = plt.gca()

        # Prepare data: Remove pairs that have zero left andd right
        x = []
        y = []
        bait_num = 0
        for i in range(0, len(en)):
            if (en[i] != 0) and (ne[i] != 0):
                x.append(en[i])
                y.append(ne[i])
                bait_num += 1

        # Get ticks
        if xy_max is None:
            xy_max = max(max(x), max(y))
        ticks, tick_labels = self.make_ticks(xy_max)

        # Scatterplot
        # -----------

        ax_s.scatter(x, y, color=color, alpha=0.5)
        ax_s.set_xticks(ticks)
        ax_s.set_xticklabels(tick_labels)
        ax_s.set_yticks(ticks)
        ax_s.set_yticklabels(tick_labels)
        ax_s.set_xlim(-xy_max / 20, xy_max + xy_max / 20)
        ax_s.set_ylim(-xy_max / 20, xy_max + xy_max / 20)

        # Histograms
        # ----------

        lim = max(x + y)
        bins = np.arange(0, lim + bin_size, bin_size)
        counts_x, bins, patches = ax_hx.hist(x, bins=bins, color=color)
        counts_y, bins, patches = ax_hy.hist(y, bins=bins, orientation='horizontal', color=color)
        xy_hist_max = max(max(counts_x), max(counts_y))

        # set x-axes
        ax_hx.set_xticks(ticks)
        ax_hy.set_yticks(ticks)
        ax_hx.set_xlim(ax_s.get_xlim())
        ax_hy.set_ylim(ax_s.get_ylim())

        # set y-axes
        ax_hx.set_yticks([0, xy_hist_max])
        ax_hy.set_xticks([0, xy_hist_max])
        ax_hx.set_ylim(0, xy_hist_max)
        ax_hy.set_xlim(0, xy_hist_max)

        ax_hx.set_title(i_cat + ' (n=' + "{:,}".format(bait_num) + ')', loc='left')
        ax_s.set_xlabel('EN')
        ax_s.set_ylabel('NE')

        return ax_hx, ax_hy, ax_s

############################


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