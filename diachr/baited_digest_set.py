from .baited_digest import BaitedDigest
from .diachromatic_interaction_set import DiachromaticInteractionSet
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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
                    print("\t[INFO] Read " + "{:,}".format(
                        self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + " interactions ...")

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
            print("\t[INFO] Total number of interactions read: " + "{:,}".format(
                self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']))
            print("\t[INFO] Total number of baited digests: " + "{:,}".format(
                self._ingest_interaction_set_info_dict['BAITED_DIGESTS']))
            print("[INFO] ... done.")

        return self._ingest_interaction_set_info_dict

    def get_ingest_interaction_set_info_report(self):
        """
        :return: Formatted string with information about ingestion of interactions for output on the screen or in report
         files.
        """

        report = "[INFO] Report on ingestion of interactions:" + '\n'
        report += "\t[INFO] Total number of interactions read: " + "{:,}".format(
            self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + '\n'
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
        report += "\t[INFO] Total number of baited digests: " + "{:,}".format(
            self._ingest_interaction_set_info_dict['BAITED_DIGESTS']) + '\n'
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

    def add_di_simple_twisted_categories_for_all_baits(self):
        """
        Add simple and twisted categories for all baits.
        """
        for chrom in self._baited_digest_dict.keys():
            for baited_digest in self._baited_digest_dict[chrom].values():
                baited_digest.add_di_simple_twisted_categories()

    # ANALYZE INTERACTION DISTANCES INDEPENDENTLY OF BAITS

    def get_median_interaction_distance(self, i_cat: str = None, e_cat: str = None):
        """
        This function calculates the median interaction distance over all baits.
        :param i_cat:
        :param e_cat:
        :return:
        """

        # Iterate chromosomes and baits on chromosomes
        i_dist_list = []
        for chrom in self._baited_digest_dict.keys():
            for key, baited_digest in self._baited_digest_dict[chrom].items():
                i_dist_list.extend(baited_digest.get_i_dist_list(i_cat, e_cat))
        return np.median(i_dist_list)

    def get_mean_interaction_distance(self, i_cat: str = None, e_cat: str = None):
        """
        This function calculates the mean interaction distance over all baits.
        :param i_cat:
        :param e_cat:
        :return:
        """

        # Iterate chromosomes and baits on chromosomes
        i_dist_list = []
        for chrom in self._baited_digest_dict.keys():
            for key, baited_digest in self._baited_digest_dict[chrom].items():
                i_dist_list.extend(baited_digest.get_i_dist_list(i_cat, e_cat))
        return np.mean(i_dist_list)

    def get_empty_pair_dict(self,
                            i_cats: list = None,
                            i_cat_colors: list = None,
                            i_cat_names: list = None):
        """
        This function initializes a data structure that can be filled with interaction numbers, read pair
        numbers or median interaction distances.

        The data structure consists of a dictionary with further sub-dictionaries. At the top level there are four
        dictionaries, one for each interaction category. In the level below, there are two lists for each category,
        one for NE and one for EN.

        In addition to the lists with the pairs, there is a list with chromosomes that were taken into account.

        When number pairs are stored in this data structure, for a given interaction category, two list elements with
        the same index should form a pair. In this case, the lists for NE and EN are the same length for each
        interaction category.

        The data structure can be filled with number pairs using the following functions 'get_*_pairs_at_baits()' and a
        filled 'pair_dict' can be plotted using the function 'get_pair_scatter_plots_with_histograms()'.

        :return: A dictionary containing lists of numbers for the various interaction aand enrichment categories and
        a list of chromosomes that were taken into account.
        """

        pair_dict = dict()

        e_cats = ['NE', 'EN']
        for i in range(len(i_cats)):
            pair_dict[i_cats[i]] = dict()
            for e_cat in e_cats:
                pair_dict[i_cats[i]][e_cat] = []
            if i_cat_names is not None:
                pair_dict[i_cats[i]]['NAME'] = i_cat_names[i]
            if i_cat_colors is not None:
                pair_dict[i_cats[i]]['COLOR'] = i_cat_colors[i]
        pair_dict['CHROMOSOMES'] = []

        return pair_dict

    def get_all_rp_nums_or_i_dists(self,
                                   chromosomes: [str] = None,
                                   number_type: str = None,
                                   verbose: bool = False):
        """
        This function collects read pair numbers or interaction distances in the various interaction and enrichment categories.


        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21']
        :param verbose: If true, messages about progress will be written to the screen
        :param number_type: Number type must be either 'RP_NUM' or 'I_DIST'
        :return: A dictionary containing the results and a list of chromosomes that were taken into account.
        """

        if (number_type is None) or number_type not in ['RP_NUM', 'I_DIST']:
            print("[ERROR] Invalid number pair type! Use one of the following:")
            print("\t[ERROR] 'RP_NUM' - Read pair numbers")
            print("\t[ERROR] 'I_DIST' - Interaction distances")
            return 1

        # Prepare data structure for results
        num_dict = self.get_empty_pair_dict(
            i_cats = ['DI', 'DI_S', 'DI_T', 'UIR', 'UI', 'ALL'],
            i_cat_names = ['Directed',
                           'Directed simple',
                           'Directed twisted',
                           'Undirected reference 1',
                           'Undirected',
                           'Undirected reference 2'],
            i_cat_colors = ['orange',
                            'pink',
                            'cadetblue',
                            'lightblue',
                            'lightgray',
                            'cornflowerblue'])
        if number_type == 'RP_NUM':
            num_dict['NUM_TYPE'] = 'Read pair number'
        if number_type == 'I_DIST':
            num_dict['NUM_TYPE'] = 'Interaction distance'

        if verbose:
            print("[INFO] Getting all " + num_dict['NUM_TYPE'].lower() + "s ...")

        # Combine lists of distances or read pair numbers from all baits
        for chrom in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chrom not in chromosomes:
                    continue

            if verbose:
                print("\t[INFO] Processing chromosome " + chrom + " ...")

            num_dict['CHROMOSOMES'].append(chrom)
            self.add_di_simple_twisted_categories_for_all_baits()
            for baited_digest in self._baited_digest_dict[chrom].values():
                for i_cat in ['DI', 'DI_S', 'DI_T', 'UIR', 'UI', 'ALL']:
                    for e_cat in ['NE', 'EN']:
                        if number_type == 'RP_NUM':
                            num_dict[i_cat][e_cat].extend(baited_digest.get_rp_num_list(i_cat, e_cat))
                        if number_type == 'I_DIST':
                            num_dict[i_cat][e_cat].extend(baited_digest.get_i_dist_list(i_cat, e_cat))

        if verbose:
            print("[INFO] ... done.")

        return num_dict

    def get_all_rp_nums_or_i_dists_histograms(self,
                                              num_dict: dict = None,
                                              y_max_set = None,
                                              q_lim: float = 0.95,
                                              description: str = "DESCRIPTION",
                                              pdf_file_name: str = "num_histograms.pdf"):
        """
        This function creates the histograms for read pair numbers or interaction distances in the various
        interaction and enrichment categories.

        :param y_max_set: Uniform maximum value for the y-axes of all histograms
        :param num_dict: Dictionary that was created with the function 'get_all_rp_nums_or_i_dists()'
        :param q_lim: Chose upper limit for x-axes based on quantile (affects presentation only)
        :param description: Brief description that is shown in the plot above the histograms
        :param pdf_file_name: Name of the PDF file that will be created.
        :return: A matplotlib 'Figure' object that can be displayed in Jupyter notebooks
        """

        # Interaction categories
        i_cats = ['DI', 'UIR', 'UI', 'ALL']

        # Prepare grid for individual plots
        fig, ax = plt.subplots(nrows=5, ncols=2, figsize=(9.5, 11),
                               gridspec_kw={'height_ratios': [0.5, 1, 1, 1, 1]})

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
        fig.text(0.015, 0.97, num_dict['NUM_TYPE'] + 's', fontsize=18,
                 fontweight='bold')
        fig.text(0.030, 0.94, 'Description: ' + description, fontsize=12)
        fig.text(0.030, 0.92, 'For chromosomes:', fontsize=12)
        if len(num_dict['CHROMOSOMES']) < 22:
            fig.text(0.045, 0.90, '[' + ", ".join(i for i in num_dict['CHROMOSOMES']) + ']', fontsize=8)
        else:
            fig.text(0.045, 0.90, '[' + ", ".join(i for i in num_dict['CHROMOSOMES'][:22]) + ',', fontsize=8)
            fig.text(0.045, 0.88, ", ".join(i for i in num_dict['CHROMOSOMES'][22:]) + ']', fontsize=8)

        # Set variables that all histograms have in common
        x_lim = 0
        for i in [0, 1, 2, 3]:
            q = max(np.quantile(num_dict[i_cats[i]]['NE'], q_lim), np.quantile(num_dict[i_cats[i]]['EN'], q_lim))
            if x_lim < q:
                x_lim = q
        x_ticks, x_tick_labels = self.make_ticks(x_lim)
        bin_width = int(x_lim / 30)

        # Create two histograms for each interaction category
        # ---------------------------------------------------

        # Prepare bins
        x_max = 0
        for i in [0, 1, 2, 3]:
            if x_max < max(num_dict[i_cats[i]]['NE'] + num_dict[i_cats[i]]['EN']):
                x_max = max(num_dict[i_cats[i]]['NE'] + num_dict[i_cats[i]]['EN'])
        bins = range(0, x_max + bin_width, bin_width)

        for i in [0, 1, 2, 3]:

            # Create histogram for NE
            counts_1, bins, patches = ax[i + 1][0].hist(
                num_dict[i_cats[i]]['NE'],
                bins=bins, density=False,
                facecolor=num_dict[i_cats[i]]['COLOR'],
                edgecolor="dimgray",
                linewidth=0.5,
                alpha=1)
            ax[i + 1][0].set_title(
                num_dict[i_cats[i]]['NAME'] + ' - NE',
                loc='left')
            ax[i + 1][0].set_xlabel(num_dict['NUM_TYPE'])
            ax[i + 1][0].set_ylabel('Frequency')
            ax[i + 1][0].set_xticks(x_ticks)
            ax[i + 1][0].set_xticklabels(x_tick_labels)
            ax[i + 1][0].set_xlim(0, x_lim)

            # Create histogram for EN
            counts_2, bins, patches = ax[i + 1][1].hist(
                num_dict[i_cats[i]]['EN'],
                bins=bins, density=False,
                facecolor=num_dict[i_cats[i]]['COLOR'],
                edgecolor="dimgray",
                linewidth=0.5,
                alpha=1)
            ax[i + 1][1].set_title(
                num_dict[i_cats[i]]['NAME'] + ' - EN',
                loc='left')
            ax[i + 1][1].set_xlabel(num_dict['NUM_TYPE'])
            ax[i + 1][1].set_ylabel('Frequency')
            ax[i + 1][1].set_xticks(x_ticks)
            ax[i + 1][1].set_xticklabels(x_tick_labels)
            ax[i + 1][1].set_xlim(0, x_lim)

            # Make y-axes comparable
            if y_max_set is None:
                y_max = max(max(counts_1), max(counts_2))
            else:
                y_max = y_max_set
            y_padding = y_max / 20
            y_max += y_padding
            y_ticks, y_tick_labels = self.make_ticks(y_max)
            ax[i + 1][0].set_yticks(y_ticks)
            ax[i + 1][0].set_yticklabels(y_tick_labels)
            ax[i + 1][0].set_ylim(0, y_max)
            ax[i + 1][1].set_yticks(y_ticks)
            ax[i + 1][1].set_yticklabels(y_tick_labels)
            ax[i + 1][1].set_ylim(0, y_max)

            # Draw vertical lines and shaded areas for median and mad
            median_ne = np.median(num_dict[i_cats[i]]['NE'])
            mad_ne = stats.median_absolute_deviation(num_dict[i_cats[i]]['NE'])
            ax[i + 1][0].axvline(median_ne, linestyle='--', linewidth=0.75, color='blue', zorder=2)
            ax[i + 1][0].axvspan(median_ne, median_ne + mad_ne, color='green', alpha=0.25, zorder=0)
            median_en = np.median(num_dict[i_cats[i]]['EN'])
            mad_en = stats.median_absolute_deviation(num_dict[i_cats[i]]['EN'])
            ax[i + 1][1].axvline(median_en, linestyle='--', linewidth=0.75, color='blue', zorder=2)
            ax[i + 1][1].axvspan(median_en, median_en + mad_en, color='green', alpha=0.25, zorder=0)

            # Add text labels with total read pair or interaction numbers, median and median absolute deviation
            ax[i + 1][0].text(x_lim - (x_lim / 3.5),
                              y_max - (y_max / 3.25),
                              'n: ' + "{:,}".format(len(num_dict[i_cats[i]]['NE'])) + '\n' + 'Mdn: ' + "{:,.0f}".format(
                                  median_ne) + '\n' + 'Mad: ' + "{:,.0f}".format(mad_ne),
                              fontsize=9,
                              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
            ax[i + 1][1].text(x_lim - (x_lim / 3.5),
                              y_max - (y_max / 3.25),
                              'n: ' + "{:,}".format(len(num_dict[i_cats[i]]['EN'])) + '\n' + 'Mdn: ' + "{:,.0f}".format(
                                  median_en) + '\n' + 'Mad: ' + "{:,.0f}".format(mad_en),
                              fontsize=9,
                              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        # Save and return figure
        fig.tight_layout(pad=1.5)
        fig.savefig(pdf_file_name)
        return fig

    def get_all_rp_nums_or_i_dists_denisty_diff_plot(self,
                                                     num_dict: dict = None,
                                                     i_cats: list = None,
                                                     q_lim: float = 0.95,
                                                     description: str = "DESCRIPTION",
                                                     pdf_file_name: str = "density_diff_plot.pdf"):
        """
        This function creates the density difference plot for a given pair of interaction categories.

        :param num_dict: A data structure that contains read pair numbers and distances for all interaction categories
        :param i_cats: List of interaction categories to be compared, e.g. ['DI', 'UIR']
        :param q_lim: Upper limit at quantile level, e.g. at 0.95, 5% of the data with the larger quantile are not shown
        :param description: Brief description that is shown in the plot above the histograms
        :param pdf_file_name:
        :return:
        """

        # Add verbose output
        # Catch wrong input, e.g. wron size of 'i_cats' list

        # Prepare grid for individual plots
        fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(14.25, 8.55),
                               gridspec_kw={'height_ratios': [0.5, 1, 1, 1]})

        # Hide unnecessary subplots
        ax[0][2].axis('off')
        ax[3][2].axis('off')

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
        fig.text(0.015*(9.5/14.25), 0.960, num_dict['NUM_TYPE'] + 's', fontsize=18,
                 fontweight='bold')
        fig.text(0.030*(9.5/14.25), 0.921, 'Description: ' + description, fontsize=12)
        fig.text(0.030*(9.5/14.25), 0.895, 'For chromosomes:', fontsize=12)
        if len(num_dict['CHROMOSOMES']) < 22:
            fig.text(0.045*(9.5/14.25), 0.869, '[' + ", ".join(i for i in num_dict['CHROMOSOMES']) + ']', fontsize=8)
        else:
            fig.text(0.045*(9.5/14.25), 0.869, '[' + ", ".join(i for i in num_dict['CHROMOSOMES'][:22]) + ',', fontsize=8)
            fig.text(0.045*(9.5/14.25), 0.843, ", ".join(i for i in num_dict['CHROMOSOMES'][22:]) + ']', fontsize=8)

        # Set variables that all histograms have in common
        x_lim = 0
        for i in [0, 1]:
            q = max(np.quantile(num_dict[i_cats[i]]['NE'], q_lim), np.quantile(num_dict[i_cats[i]]['EN'], q_lim))
            if x_lim < q:
                x_lim = q
        x_ticks, x_tick_labels = self.make_ticks(x_lim)
        bin_width = int(x_lim / 30)

        # Prepare bins
        x_max = 0
        for i in [0, 1]:
            if x_max < max(num_dict[i_cats[i]]['NE'] + num_dict[i_cats[i]]['EN']):
                x_max = max(num_dict[i_cats[i]]['NE'] + num_dict[i_cats[i]]['EN'])
        bins = range(0, x_max + bin_width, bin_width)

        # Create histograms for the two categories
        # ----------------------------------------

        abs_ne = []
        abs_en = []
        densities_ne = []
        densities_en = []
        for i in [0, 1]:

            # Create histogram for NE
            counts_ne, bins, patches = ax[i + 1][0].hist(
                num_dict[i_cats[i]]['NE'],
                bins=bins, density=False,
                facecolor=num_dict[i_cats[i]]['COLOR'],
                edgecolor="dimgray",
                linewidth=0.5,
                alpha=1)
            ax[i + 1][0].set_title(num_dict[i_cats[i]]['NAME'] + ' - NE', loc='left')
            ax[i + 1][0].set_xlabel(num_dict['NUM_TYPE'])
            ax[i + 1][0].set_ylabel('Frequency')
            ax[i + 1][0].set_xticks(x_ticks)
            ax[i + 1][0].set_xticklabels(x_tick_labels)
            ax[i + 1][0].set_xlim(0, x_lim)

            # Create histogram for NE
            counts_en, bins, patches = ax[i + 1][1].hist(
                num_dict[i_cats[i]]['EN'],
                bins=bins, density=False,
                facecolor=num_dict[i_cats[i]]['COLOR'],
                edgecolor="dimgray",
                linewidth=0.5,
                alpha=1)
            ax[i + 1][1].set_title(num_dict[i_cats[i]]['NAME'] + ' - EN', loc='left')
            ax[i + 1][1].set_xlabel(num_dict['NUM_TYPE'])
            ax[i + 1][1].set_ylabel('Frequency')
            ax[i + 1][1].set_xticks(x_ticks)
            ax[i + 1][1].set_xticklabels(x_tick_labels)
            ax[i + 1][1].set_xlim(0, x_lim)

            # Draw vertical lines and shaded areas for median and MAD
            median_ne = np.median(num_dict[i_cats[i]]['NE'])
            mad_ne = stats.median_absolute_deviation(num_dict[i_cats[i]]['NE'])
            ax[i + 1][0].axvline(median_ne, linestyle='--', linewidth=0.75, color='blue', zorder=2)
            ax[i + 1][0].axvspan(median_ne, median_ne + mad_ne, color='green', alpha=0.25, zorder=0)
            median_en = np.median(num_dict[i_cats[i]]['EN'])
            mad_en = stats.median_absolute_deviation(num_dict[i_cats[i]]['EN'])
            ax[i + 1][1].axvline(median_en, linestyle='--', linewidth=0.75, color='blue', zorder=2)
            ax[i + 1][1].axvspan(median_en, median_en + mad_en, color='green', alpha=0.25, zorder=0)

            # Keep track of bin counts and denisties densities
            abs_ne.append(counts_ne)
            abs_en.append(counts_en)
            densities_ne.append([count / sum(counts_ne) for count in counts_ne])
            densities_en.append([count / sum(counts_en) for count in counts_en])

        # Transform bins to center bin positions. The original list for the bins contains the boundaries between bins,
        # but for the bar plot we need the center positions of bins.
        bcp = bin_width / 2
        bcp_list = [bcp]
        for i in range(len(densities_ne[0]) - 1):
            bcp += bin_width
            bcp_list.append(bcp)

        # Add second axes with densities and normalize all histograms to maximum density
        # -----------------------------------------------------------------------------

        # Determine maximal density in all four plots
        yd_max = max(max(densities_ne[0]), max(densities_en[0]), max(densities_ne[1]), max(densities_en[1]))
        y_padding = yd_max / 20
        yd_max += y_padding

        # Add density axes, normalize and add text labels
        for i in [0, 1]:

            # Left histogram (NNE)
            yc_max = yd_max * sum(abs_ne[i])
            yc_ticks, yc_tick_labels = self.make_ticks(yc_max)
            ax[i + 1][0].set_yticks(yc_ticks)
            ax[i + 1][0].set_yticklabels(yc_tick_labels)
            ax[i + 1][0].set_ylim(0, yc_max)
            ax_dens = ax[i + 1][0].twinx()
            ax_dens.plot(bcp_list, densities_ne[i], color='red', linewidth=0.5)
            ax_dens.set_ylabel('Density')
            ax_dens.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
            ax_dens.set_ylim(0, yd_max)
            # Add text labels with total read pair or interaction numbers, median and median absolute deviation
            ax[i + 1][0].text(x_lim - (x_lim / 3),
                              yc_max - (yc_max / 3.2),
                              'n: ' + "{:,}".format(len(num_dict[i_cats[i]]['NE'])) + '\n' + 'Mdn: ' + "{:,.0f}".format(
                                  np.median(num_dict[i_cats[i]]['NE'])) + '\n' + 'Mad: ' + "{:,.0f}".format(stats.median_absolute_deviation(num_dict[i_cats[i]]['NE'])),
                              fontsize=9,
                              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

            # Right histogram (EN)
            yc_max = yd_max * sum(abs_en[i])
            yc_ticks, yc_tick_labels = self.make_ticks(yc_max)
            ax[i + 1][1].set_yticks(yc_ticks)
            ax[i + 1][1].set_yticklabels(yc_tick_labels)
            ax[i + 1][1].set_ylim(0, yc_max)
            ax_dens = ax[i + 1][1].twinx()
            ax_dens.plot(bcp_list, densities_en[i], color='red', linewidth=0.5)
            ax_dens.set_ylabel('Density')
            ax_dens.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
            ax_dens.set_ylim(0, yd_max)
            # Add text labels with total read pair or interaction numbers, median and median absolute deviation
            ax[i + 1][1].text(x_lim - (x_lim / 3),
                              yc_max - (yc_max / 3.2),
                              'n: ' + "{:,}".format(len(num_dict[i_cats[i]]['EN'])) + '\n' +
                              'Mdn: ' + "{:,.0f}".format(np.median(num_dict[i_cats[i]]['EN'])) + '\n' +
                              'Mad: ' + "{:,.0f}".format(stats.median_absolute_deviation(num_dict[i_cats[i]]['EN'])),
                              fontsize=9,
                              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        # Add subplots for density differences
        # ------------------------------------

        # Get four lists with densities differences from all four density lists
        density_diff_cat_0_cat_1_ne = [] # Below histograms on the left
        density_diff_cat_0_cat_1_en = [] # Below histograms on the right
        density_diff_cat_0_ne_en = []    # To the left of the upper histograms
        density_diff_cat_1_ne_en = []    # To the left of the lower histograms
        for i in range(0, len(densities_ne[0])):
            density_diff_cat_0_cat_1_ne.append(densities_ne[0][i] - densities_ne[1][i])
            density_diff_cat_0_cat_1_en.append(densities_en[0][i] - densities_en[1][i])
            density_diff_cat_0_ne_en.append(densities_ne[0][i] - densities_en[0][i])
            density_diff_cat_1_ne_en.append(densities_ne[1][i] - densities_en[1][i])

        # Determine sum of density differences
        dd_cat_0_cat_1_ne_sum = sum(map(abs, density_diff_cat_0_cat_1_ne))
        dd_cat_0_cat_1_en_sum = sum(map(abs, density_diff_cat_0_cat_1_en))
        dd_cat_0_ne_en_sum = sum(map(abs, density_diff_cat_0_ne_en))
        dd_cat_1_ne_en_sum = sum(map(abs, density_diff_cat_1_ne_en))

        # Bar plot below histograms on the left
        ax[3][0].bar(bcp_list,
                     density_diff_cat_0_cat_1_ne,
                     width=bin_width,
                     color=[num_dict[i_cats[0]]['COLOR'] if 0 < dd else num_dict[i_cats[1]]['COLOR'] for dd in density_diff_cat_0_cat_1_ne],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[3][0].set_title('NE', loc='left')
        ax[3][0].set_xticks(x_ticks)
        ax[3][0].set_xticklabels(x_tick_labels)
        ax[3][0].set_xlim(0, x_lim)
        ax[3][0].set_xlabel(num_dict['NUM_TYPE'])
        ax[3][0].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[3][0].ticklabel_format(axis='y', style='sci', scilimits=(-2,0), useMathText=False)
        ax[3][0].yaxis.tick_right()
        ax[3][0].yaxis.set_ticks_position('both')
        ax[3][0].set_ylabel('Density difference')
        ax[3][0].yaxis.set_label_position('right')

        # Bar plot below histograms on the right
        ax[3][1].bar(bcp_list,
                     density_diff_cat_0_cat_1_en,
                     width=bin_width,
                     color=[num_dict[i_cats[0]]['COLOR'] if 0 < dd else num_dict[i_cats[1]]['COLOR'] for dd in density_diff_cat_0_cat_1_en],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[3][1].set_title('EN', loc='left')
        ax[3][1].set_xticks(x_ticks)
        ax[3][1].set_xticklabels(x_tick_labels)
        ax[3][1].set_xlim(0, x_lim)
        ax[3][1].set_xlabel(num_dict['NUM_TYPE'])
        ax[3][1].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[3][1].ticklabel_format(axis='y', style='sci', scilimits=(-2,0))
        ax[3][1].yaxis.tick_right()
        ax[3][1].yaxis.set_ticks_position('both')
        ax[3][1].set_ylabel('Density difference')
        ax[3][1].yaxis.set_label_position('right')

        # Make y-axes comparable for the two bar plots below the histograms and add sums of density differences
        y_min = min(min(density_diff_cat_0_cat_1_ne), min(density_diff_cat_0_cat_1_en))
        y_max = max(max(density_diff_cat_0_cat_1_ne), max(density_diff_cat_0_cat_1_en))
        y_padding = ((y_max - y_min)/20)
        y_min -= y_padding
        y_max += y_padding
        ax[3][0].set_ylim(y_min, y_max)
        ax[3][1].set_ylim(y_min, y_max)
        ax[3][0].text(x_lim - (x_lim / 3.2),
                      y_max - ((y_max - y_min) / 8),
                      'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_cat_1_ne_sum),
                      fontsize=9,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
        ax[3][1].text(x_lim - (x_lim / 3.2),
                      y_max - ((y_max - y_min) / 8),
                      'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_cat_1_en_sum),
                      fontsize=9,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        # Bar plot to the left of the upper histograms
        ax[1][2].bar(bcp_list,
                     density_diff_cat_0_ne_en,
                     width=bin_width,
                     color=[num_dict[i_cats[0]]['COLOR'] if 0 < dd else num_dict[i_cats[0]]['COLOR'] for dd in density_diff_cat_0_ne_en],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[1][2].set_title( num_dict[i_cats[0]]['NAME'], loc='left')
        ax[1][2].set_xticks(x_ticks)
        ax[1][2].set_xticklabels(x_tick_labels)
        ax[1][2].set_xlim(0, x_lim)
        ax[1][2].set_xlabel(num_dict['NUM_TYPE'])
        ax[1][2].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[1][2].ticklabel_format(axis='y', style='sci', scilimits=(-2, 0), useMathText=False)
        ax[1][2].yaxis.tick_right()
        ax[1][2].yaxis.set_ticks_position('both')
        ax[1][2].set_ylabel('Density difference')
        ax[1][2].yaxis.set_label_position('right')

        # Bar plot to the left of the lower histograms
        ax[2][2].bar(bcp_list,
                     density_diff_cat_1_ne_en,
                     width=bin_width,
                     color=[num_dict[i_cats[1]]['COLOR'] if 0 < dd else num_dict[i_cats[1]]['COLOR'] for dd in density_diff_cat_1_ne_en],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[2][2].set_title(num_dict[i_cats[1]]['NAME'], loc='left')
        ax[2][2].set_xticks(x_ticks)
        ax[2][2].set_xticklabels(x_tick_labels)
        ax[2][2].set_xlim(0, x_lim)
        ax[2][2].set_xlabel(num_dict['NUM_TYPE'])
        ax[2][2].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[2][2].ticklabel_format(axis='y', style='sci', scilimits=(-2, 0), useMathText=False)
        ax[2][2].yaxis.tick_right()
        ax[2][2].yaxis.set_ticks_position('both')
        ax[2][2].set_ylabel('Density difference')
        ax[2][2].yaxis.set_label_position('right')

        # Make y-axes comparable for the two bar plots to the left of the histograms and add sums of density differences
        y_min = min(min(density_diff_cat_0_ne_en), min(density_diff_cat_1_ne_en))
        y_max = max(max(density_diff_cat_0_ne_en), max(density_diff_cat_1_ne_en))
        y_padding = ((y_max - y_min) / 20)
        y_min -= y_padding
        y_max += y_padding
        ax[1][2].set_ylim(y_min, y_max)
        ax[2][2].set_ylim(y_min, y_max)
        ax[1][2].text(x_lim - (x_lim / 3.2),
                      y_max - ((y_max - y_min) / 8),
                      'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_ne_en_sum),
                      fontsize=9,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
        ax[2][2].text(x_lim - (x_lim / 3.2),
                      y_max - ((y_max - y_min) / 8),
                      'sum(|dd|): ' + "{:.2f}".format(dd_cat_1_ne_en_sum),
                      fontsize=9,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        # Save and return figure
        fig.tight_layout(pad=1.5)
        fig.savefig(pdf_file_name)
        return fig

    # PAIRS OF NUMBERS AT BAITS

    def get_number_pairs_at_baits(self,
                                  number_pair_type: str = None,
                                  chromosomes: [str] = None,
                                  verbose: bool = False):
        """
        This function determines a pair of numbers for each bait, the number to the left and the number to the right.
        Numbers can be determined for interaction numbers, read pair numbers and median interaction distances. The
        results are returned in a special data structure called 'pair_dict'.

        :param number_pair_type: Type of distance pairs that will be determined either 'I_NUM', 'RP_NUM' or 'MED_I_DIST'.
        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21'].
        :param verbose: If true, messages about progress will be written to the screen.
        :return: A dictionary containing the number pairs (see 'get_empty_pair_dict()').
        """

        if (number_pair_type is None) or number_pair_type not in ['I_NUM', 'RP_NUM', 'MED_RP_NUM', 'MED_I_DIST']:
            print("[ERROR] Invalid number pair type! Use one of the following:")
            print("\t[ERROR] 'I_NUM' - Interaction numbers")
            print("\t[ERROR] 'RP_NUM' - Read pair numbers")
            print("\t[ERROR] 'MED_RP_NUM' - Median read pair numbers")
            print("\t[ERROR] 'MED_I_DIST' - Median interaction distances")
            return 1

        # Prepare data structure for results
        i_cats = ['DI', 'UIR', 'UI', 'ALL']
        pair_dict = self.get_empty_pair_dict(['DI', 'UIR', 'UI', 'ALL'])
        pair_dict['BAIT_NUM_TOTAL'] = 0
        if number_pair_type == 'I_NUM':
            pair_dict['NUM_PAIR_TYPE'] = 'Interaction number'
        if number_pair_type == 'RP_NUM':
            pair_dict['NUM_PAIR_TYPE'] = 'Read pair number'
        if number_pair_type == 'MED_RP_NUM':
            pair_dict['NUM_PAIR_TYPE'] = 'Median read pair number'
        if number_pair_type == 'MED_I_DIST':
            pair_dict['NUM_PAIR_TYPE'] = 'Median interaction distance'

        if verbose:
            print("[INFO] Determining pairs of " + pair_dict[
                'NUM_PAIR_TYPE'].lower() + "s (to the left and right) at baits ...")

        # Iterate chromosomes and baits on chromosomes
        for chrom in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chrom not in chromosomes:
                    continue

            if verbose:
                print("\t[INFO] Processing chromosome " + chrom + " ...")

            # Determine number pairs for individual baits
            pair_dict['CHROMOSOMES'].append(chrom)
            for key, baited_digest in self._baited_digest_dict[chrom].items():
                pair_dict['BAIT_NUM_TOTAL'] += 1
                for i_cat in i_cats:
                    if number_pair_type == 'I_NUM':
                        pair_dict[i_cat]['NE'].append(baited_digest.get_interaction_number(i_cat, 'NE'))
                        pair_dict[i_cat]['EN'].append(baited_digest.get_interaction_number(i_cat, 'EN'))
                    if number_pair_type == 'RP_NUM':
                        pair_dict[i_cat]['NE'].append(baited_digest.get_read_pair_number(i_cat, 'NE'))
                        pair_dict[i_cat]['EN'].append(baited_digest.get_read_pair_number(i_cat, 'EN'))
                    if number_pair_type == 'MED_RP_NUM':
                        pair_dict[i_cat]['NE'].append(baited_digest.get_median_read_pair_number(i_cat, 'NE'))
                        pair_dict[i_cat]['EN'].append(baited_digest.get_median_read_pair_number(i_cat, 'EN'))
                    if number_pair_type == 'MED_I_DIST':
                        pair_dict[i_cat]['NE'].append(baited_digest.get_median_interaction_distance(i_cat, 'NE'))
                        pair_dict[i_cat]['EN'].append(baited_digest.get_median_interaction_distance(i_cat, 'EN'))

        if verbose:
            print("[INFO] ... done.")

        return pair_dict

    def make_ticks(self, max_val: int = 100):

        if max_val < 10:
            return list(range(0, 50 + 1, 1)), list(range(0, 50 + 1, 1))

        if max_val < 25:
            return list(range(0, 50 + 1, 5)), list(range(0, 50 + 1, 5))

        if max_val < 50:
            return list(range(0, 50 + 1, 10)), list(range(0, 50 + 1, 10))

        # Create list of maximal ticks
        tick_lims = []
        x = [10, 25, 50]
        while x[2] <= max_val:
            x = [(y * 10) for y in x]
            tick_lims = tick_lims + x

        # Find maximal tick for input value
        tick_max_idx = 0
        while tick_lims[tick_max_idx] <= max_val:
            tick_max_idx += 1

        ticks = list(range(0, int(tick_lims[tick_max_idx]) + 1, int(tick_lims[tick_max_idx] / 5)))

        if max_val < 1000:
            tick_labels = ticks
        elif max_val < 1000000:
            tick_labels = []
            for tick in ticks:
                tick_labels.append(str(tick / 1000) + 'k')
        else:
            tick_labels = []
            for tick in ticks:
                tick_labels.append(str(tick / 1000000) + 'M')

        return ticks, tick_labels

    def create_single_pair_scatter_plot_with_histograms(self,
                                                        i_cat_label: str = 'I_CAT_LABEL',
                                                        i_cat_color: str = 'black',
                                                        en_list: list() = None,
                                                        ne_list: list() = None,
                                                        x_lab: str = 'EN',
                                                        y_lab: str = 'NE',
                                                        bin_size: int = 1,
                                                        xy_max: int = None,
                                                        add_text_labels=True,
                                                        draw_mean_and_sd=True,
                                                        ax_hx=None,
                                                        ax_hy=None,
                                                        ax_s=None,
                                                        **plt_kwargs):
        """
        This function generates a scatter plot, which is supplemented by histograms on the x and y axes.

        :param i_cat_label: Label for interaction category
        :param i_cat_color: Color for interaction category
        :param en_list: List of numbers that pair with the numbers in the other list of the same size
        :param ne_list: Other list of numbers (en_list[i] and ne_list[i] form a pair)
        :param x_lab: Label for the x-axis
        :param y_lab: Label for the y-axis
        :param bin_size: Bin size for histograms on axes
        :param xy_max: Maximum value on x and y axes
        :param add_text_labels: If true, then add text labels with bait and interaction counts
        :param draw_mean_and_sd: If true, dashed line and shaded area for mean and standard deviations will be added
        :param ax_hx: 'matplotlib' object for histogram along the x axis
        :param ax_hy: 'matplotlib' object for histogram along the y axis
        :param ax_s: 'matplotlib' object for scatter plot
        :param plt_kwargs: Required for this function to work
        :return: The three 'matplotlib' objects filled with content
        """

        # Required for this function to work
        if ax_hx is None:
            ax_hx = plt.gca()
        if ax_hy is None:
            ax_hy = plt.gca()
        if ax_s is None:
            ax_s = plt.gca()

        # Prepare data: Remove pairs that have zero left and right
        x = []
        y = []
        bait_num = 0
        i_num_en_total = 0
        i_num_ne_total = 0
        for i in range(0, len(en_list)):
            if not (en_list[i] == 0 and ne_list[i] == 0):
                x.append(en_list[i])
                y.append(ne_list[i])
                i_num_en_total += en_list[i]
                i_num_ne_total += ne_list[i]
                bait_num += 1

        # Get ticks
        if xy_max is None:
            xy_max = max(max(x), max(y))
        ticks, tick_labels = self.make_ticks(xy_max)

        # Scatter plot
        # ------------

        ax_s.scatter(x, y, color=i_cat_color, alpha=0.5)
        ax_s.set_xticks(ticks)
        ax_s.set_xticklabels(tick_labels)
        ax_s.set_yticks(ticks)
        ax_s.set_yticklabels(tick_labels)
        ax_s.set_xlim(-xy_max / 20, xy_max + xy_max / 20)
        ax_s.set_ylim(-xy_max / 20, xy_max + xy_max / 20)

        # Add line and area for mean and standard deviation
        if draw_mean_and_sd:
            mean_en = np.mean(en_list)
            std_en = np.std(en_list)
            ax_s.axvspan(mean_en, mean_en + std_en, color='green', alpha=0.25, zorder=0)
            ax_s.axvline(mean_en, linestyle='--', linewidth=0.5, color='blue')
            mean_ne = np.mean(ne_list)
            std_ne = np.std(ne_list)
            ax_s.axhspan(mean_ne, mean_ne + std_ne, color='green', alpha=0.25, zorder=0)
            ax_s.axhline(mean_ne, linestyle='--', linewidth=0.5, color='blue')

        # Histograms
        # ----------

        lim = max(x + y)
        bins = np.arange(0, lim + bin_size, bin_size)
        counts_x, bins, patches = ax_hx.hist(x, bins=bins, color=i_cat_color)
        counts_y, bins, patches = ax_hy.hist(y, bins=bins, orientation='horizontal', color=i_cat_color)
        xy_hist_max = max(max(counts_x), max(counts_y))

        # set x-axes
        ax_hx.ticklabel_format(useOffset=False, style='plain')
        ax_hx.set_xticks(ticks)
        ax_hy.set_yticks(ticks)
        ax_hx.set_xlim(ax_s.get_xlim())
        ax_hy.set_ylim(ax_s.get_ylim())

        # set y-axes
        ax_hy.ticklabel_format(useOffset=False, style='plain')
        ax_hx.set_yticks([0, xy_hist_max])
        ax_hy.set_xticks([0, xy_hist_max])
        ax_hx.set_ylim(0, xy_hist_max)
        ax_hy.set_xlim(0, xy_hist_max)

        # Add text labels
        if add_text_labels:
            ax_s.text(xy_max - xy_max / 3.5,
                      xy_max - (xy_max / 15.5),
                      'n: ' + "{:,}".format(bait_num),
                      bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
            ax_hx.text(xy_max / 5,
                       xy_hist_max - (xy_hist_max / 4),
                       'EN: ' + "{:,}".format(i_num_en_total),
                       bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

            ax_hy.text(xy_hist_max - (xy_hist_max / 3.2),
                       xy_max / 4.4,
                       'NE: ' + "{:,}".format(i_num_ne_total),
                       rotation=-90,
                       bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        ax_hx.set_title(i_cat_label, size='x-large', loc='left')
        ax_s.set_xlabel(x_lab)
        ax_s.set_ylabel(y_lab)

        return ax_hx, ax_hy, ax_s

    def get_pair_scatter_plots_with_histograms(self,
                                               pairs_dict=None,
                                               set_xy_max=None,
                                               draw_mean_and_sd=False,
                                               sup_title: str = 'SUP_TITLE',
                                               description: str = 'DESCRIPTION',
                                               pdf_file_name: str = 'pair_scatter_plots_with_histograms.pdf'):
        """
        This function generates a figure with a text header and four scatter plots supplemented by histograms on the x
        and y axes.

        :param pairs_dict:
        :param set_xy_max:
        :param sup_title:
        :param description:
        :param pdf_file_name:
        :return:
        """

        # Define layout of the entire figure (grid of rectangles)
        # -------------------------------------------------------

        header_height = 0.15
        fig = plt.figure(figsize=(10, (1.0 + header_height) * 10))

        # 2*a + 2*b + 2*sp2 + 3*sp1 = 1.0
        a = 0.260
        b = 0.075
        sp1 = 0.10
        sp2 = 0.015

        x1 = sp1
        x2 = x1 + a + sp2
        x3 = x2 + b + sp1
        x4 = x3 + a + sp2

        y1 = x1 * (1.0 - header_height)
        y2 = x2 * (1.0 - header_height)
        y3 = x3 * (1.0 - header_height)
        y4 = x4 * (1.0 - header_height)

        r0 = [x1, y4 + (1.0 - header_height) * b + (1.0 - header_height) * sp1, x4 - x1 + b, header_height]

        r1 = [x1, y1, a, a * (1.0 - header_height)]
        r2 = [x2, y1, b, a * (1.0 - header_height)]
        r3 = [x3, y1, a, a * (1.0 - header_height)]
        r4 = [x4, y1, b, a * (1.0 - header_height)]

        r5 = [x1, y2, a, b * (1.0 - header_height)]
        r6 = [x3, y2, a, b * (1.0 - header_height)]

        r7 = [x1, y3, a, a * (1.0 - header_height)]
        r8 = [x2, y3, b, a * (1.0 - header_height)]
        r9 = [x3, y3, a, a * (1.0 - header_height)]
        r10 = [x4, y3, b, a * (1.0 - header_height)]

        r11 = [x1, y4, a, b * (1.0 - header_height)]
        r12 = [x3, y4, a, b * (1.0 - header_height)]

        # Fill rectangular areas with content
        # -----------------------------------

        # Adjust the bin size dynamically to DI and UIR
        DI_UIR_MAX = max(max(pairs_dict['DI']['NE']),
                         max(pairs_dict['DI']['EN']),
                         max(pairs_dict['UIR']['NE']),
                         max(pairs_dict['UIR']['EN']))
        BIN_SIZE = int(np.ceil(DI_UIR_MAX / 100))

        # Get total number of baits and list of chromosomes for header
        BAITS_TOTAL = pairs_dict['BAIT_NUM_TOTAL']
        CHROMOSOMES = pairs_dict['CHROMOSOMES']

        # Create header with text
        ax0_header = plt.axes(r0)
        ax0_header.spines['left'].set_color('white')
        ax0_header.spines['right'].set_color('white')
        ax0_header.spines['top'].set_color('white')
        ax0_header.spines['bottom'].set_color('white')
        ax0_header.tick_params(axis='x', colors='white')
        ax0_header.tick_params(axis='y', colors='white')
        fig.text(0.015, 0.97, sup_title, fontsize=18, fontweight='bold')
        fig.text(0.030, 0.94, 'Description: ' + description, fontsize=12)
        fig.text(0.030, 0.92, 'Total number of baits: ' + str(BAITS_TOTAL), fontsize=12)
        fig.text(0.030, 0.90, 'Bin size: ' + str(BIN_SIZE), fontsize=12)
        fig.text(0.030, 0.88, 'For chromosomes:', fontsize=12)
        if len(CHROMOSOMES) < 22:
            fig.text(0.045, 0.86, '[' + ", ".join(i for i in CHROMOSOMES) + ']', fontsize=8)
        else:
            fig.text(0.045, 0.86, '[' + ", ".join(i for i in CHROMOSOMES[:22]) + ',', fontsize=8)
            fig.text(0.045, 0.84, ", ".join(i for i in CHROMOSOMES[22:]) + ']', fontsize=8)

        # Directed interactions (DI)
        ax1_hx = plt.axes(r11)
        ax1_hx.tick_params(direction='in', labelbottom=False)
        ax1_hy = plt.axes(r8)
        ax1_hy.tick_params(direction='in', labelleft=False)
        ax1_s = plt.axes(r7)
        ax1_s.tick_params(direction='in', top=True, right=True)
        self.create_single_pair_scatter_plot_with_histograms(
            i_cat_label='DI',
            en_list=pairs_dict['DI']['EN'],
            ne_list=pairs_dict['DI']['NE'],
            x_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - EN',
            y_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - NE',
            i_cat_color='orange',
            bin_size=BIN_SIZE,
            add_text_labels=False,
            xy_max=set_xy_max,
            draw_mean_and_sd=draw_mean_and_sd,
            ax_hx=ax1_hx,
            ax_hy=ax1_hy,
            ax_s=ax1_s)

        # Undirected reference interactions (UIR)
        ax2_hx = plt.axes(r12)
        ax2_hx.tick_params(direction='in', labelbottom=False)
        ax2_hy = plt.axes(r10)
        ax2_hy.tick_params(direction='in', labelleft=False)
        ax2_s = plt.axes(r9)
        ax2_s.tick_params(direction='in', top=True, right=True)
        self.create_single_pair_scatter_plot_with_histograms(
            i_cat_label='UIR',
            en_list=pairs_dict['UIR']['EN'],
            ne_list=pairs_dict['UIR']['NE'],
            x_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - EN',
            y_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - NE',
            i_cat_color='lightblue',
            bin_size=BIN_SIZE,
            add_text_labels=False,
            xy_max=set_xy_max,
            draw_mean_and_sd=draw_mean_and_sd,
            ax_hx=ax2_hx,
            ax_hy=ax2_hy,
            ax_s=ax2_s)

        # Make axes for DI and UIR comparable
        xy_max = int(max(ax1_s.get_xlim()[1], ax2_s.get_xlim()[1]))
        ticks, tick_labels = self.make_ticks(xy_max)
        # Scatterplots
        ax1_s.set_xticks(ticks)
        ax1_s.set_xticklabels(tick_labels)
        ax1_s.set_yticks(ticks)
        ax1_s.set_yticklabels(tick_labels)
        ax1_s.set_xlim(-xy_max / 20, xy_max + xy_max / 20)
        ax1_s.set_ylim(-xy_max / 20, xy_max + xy_max / 20)
        ax2_s.set_xticks(ticks)
        ax2_s.set_xticklabels(tick_labels)
        ax2_s.set_yticks(ticks)
        ax2_s.set_yticklabels(tick_labels)
        ax2_s.set_xlim(-xy_max / 20, xy_max + xy_max / 20)
        ax2_s.set_ylim(-xy_max / 20, xy_max + xy_max / 20)
        # Histograms
        ax1_hx.set_xticks(ticks)
        ax1_hy.set_yticks(ticks)
        ax1_hx.set_xlim(ax1_s.get_xlim())
        ax1_hy.set_ylim(ax1_s.get_ylim())

        # Add text labels
        bait_num = 0
        en_list = pairs_dict['DI']['EN']
        ne_list = pairs_dict['DI']['NE']
        i_num_en_total = 0
        i_num_ne_total = 0
        for i in range(0, len(en_list)):
            if not (en_list[i] == 0 and ne_list[i] == 0):
                i_num_en_total += en_list[i]
                i_num_ne_total += ne_list[i]
                bait_num += 1

        ax1_s.text(xy_max - xy_max / 3.5,
                   xy_max - (xy_max / 15.5),
                   'n: ' + "{:,}".format(bait_num),
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        ax1_hx.text(xy_max / 5,
                    ax1_hx.get_ylim()[1] - (ax1_hx.get_ylim()[1] / 4),
                    'EN: ' + "{:,}".format(i_num_en_total),
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        ax1_hy.text(ax1_hy.get_xlim()[1] - (ax1_hy.get_xlim()[1] / 3.2),
                    xy_max / 4.4,
                    'NE: ' + "{:,}".format(i_num_ne_total),
                    rotation=-90,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        ax2_hx.set_xticks(ticks)
        ax2_hy.set_yticks(ticks)
        ax2_hx.set_xlim(ax1_s.get_xlim())
        ax2_hy.set_ylim(ax1_s.get_ylim())

        bait_num = 0
        en_list = pairs_dict['UIR']['EN']
        ne_list = pairs_dict['UIR']['NE']
        i_num_en_total = 0
        i_num_ne_total = 0
        for i in range(0, len(en_list)):
            if not (en_list[i] == 0 and ne_list[i] == 0):
                i_num_en_total += en_list[i]
                i_num_ne_total += ne_list[i]
                bait_num += 1

        ax2_s.text(xy_max - xy_max / 3.5,
                   xy_max - (xy_max / 15.5),
                   'n: ' + "{:,}".format(bait_num),
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        ax2_hx.text(xy_max / 5,
                    ax2_hx.get_ylim()[1] - (ax2_hx.get_ylim()[1] / 4),
                    'EN: ' + "{:,}".format(i_num_en_total),
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        ax2_hy.text(ax2_hy.get_xlim()[1] - (ax2_hy.get_xlim()[1] / 3.2),
                    xy_max / 4.4,
                    'NE: ' + "{:,}".format(i_num_ne_total),
                    rotation=-90,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, boxstyle='round'))

        # Undirected interactions (UI)
        ax3_hx = plt.axes(r5)
        ax3_hx.tick_params(direction='in', labelbottom=False)
        ax3_hy = plt.axes(r2)
        ax3_hy.tick_params(direction='in', labelleft=False)
        ax3_s = plt.axes(r1)
        ax3_s.tick_params(direction='in', top=True, right=True)
        self.create_single_pair_scatter_plot_with_histograms(
            i_cat_label='UI',
            en_list=pairs_dict['UI']['EN'],
            ne_list=pairs_dict['UI']['NE'],
            x_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - EN',
            y_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - NE',
            i_cat_color='gray',
            bin_size=BIN_SIZE,
            xy_max=set_xy_max,
            draw_mean_and_sd=draw_mean_and_sd,
            ax_hx=ax3_hx,
            ax_hy=ax3_hy,
            ax_s=ax3_s)

        # All interactions (ALL)
        ax4_hx = plt.axes(r6)
        ax4_hx.tick_params(direction='in', labelbottom=False)
        ax4_hy = plt.axes(r4)
        ax4_hy.tick_params(direction='in', labelleft=False)
        ax4_s = plt.axes(r3)
        ax4_s.tick_params(direction='in', top=True, right=True)
        self.create_single_pair_scatter_plot_with_histograms(
            i_cat_label='ALL',
            en_list=pairs_dict['ALL']['EN'],
            ne_list=pairs_dict['ALL']['NE'],
            x_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - EN',
            y_lab=pairs_dict['NUM_PAIR_TYPE'] + ' - NE',
            i_cat_color='pink',
            bin_size=BIN_SIZE,
            xy_max=set_xy_max,
            draw_mean_and_sd=draw_mean_and_sd,
            ax_hx=ax4_hx,
            ax_hy=ax4_hy,
            ax_s=ax4_s)

        # Save and return figure
        fig.savefig(pdf_file_name)
        return fig

    ############################

    def select_undirected_reference_interactions_at_baits(self, lower_q: float = 0.25, upper_q: float = 0.75):

        # Iterate chromosomes and baits on chromosomes
        sorted_baited_digest_keys = self.get_baited_digest_keys_sorted_by_sta_pos()
        for chrom in self._baited_digest_dict.keys():
            # Select reference interactions at individual baits
            for baited_digest_key in sorted_baited_digest_keys[chrom]:
                baited_digest = self._baited_digest_dict[chrom][baited_digest_key]
                baited_digest.select_undirected_reference_interactions(lower_q, upper_q)

    def write_bed_files_with_baited_interactions(self, out_prefix: str = 'OUT_PREFIX', chromosomes: [str] = None):

        # Open streams for BED files
        bed_file_streams = {
            'DI': [],
            'UIR': [],
            'UI': [],
            'ALL': []
        }
        for i in range(0, 10):
            bed_stream = open(out_prefix + '_baited_interactions_' + str(i + 1) + '_di.bed', 'wt')
            bed_stream.write("track type=bed name=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - DI" + "\" description=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - DI" + "\" visibility=4 itemRgb=\"On\"" + '\n')
            bed_file_streams['DI'].append(bed_stream)
            bed_stream = open(out_prefix + '_baited_interactions_' + str(i + 1) + '_uir.bed', 'wt')
            bed_stream.write("track type=bed name=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - UIR" + "\" description=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - UIR" + "\" visibility=4 itemRgb=\"On\"" + '\n')
            bed_file_streams['UIR'].append(bed_stream)
            bed_stream = open(out_prefix + '_baited_interactions_' + str(i + 1) + '_ui.bed', 'wt')
            bed_stream.write("track type=bed name=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - UI" + "\" description=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - UI" + "\" visibility=4 itemRgb=\"On\"" + '\n')
            bed_file_streams['UI'].append(bed_stream)
            bed_stream = open(out_prefix + '_baited_interactions_' + str(i + 1) + '_all.bed', 'wt')
            bed_stream.write("track type=bed name=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - ALL" + "\" description=\"" + out_prefix + " Baited interactions " + str(i + 1) +
                             " - ALL" + "\" visibility=4 itemRgb=\"On\"" + '\n')
            bed_file_streams['ALL'].append(bed_stream)

        # Iterate chromosomes and BaitedDigest objects sorted by starting position
        sorted_baited_digest_keys = self.get_baited_digest_keys_sorted_by_sta_pos()
        # Init number for BED files with interaction at baits
        bed_file_number = 0
        for chrom in self._baited_digest_dict.keys():
            if chromosomes is not None:
                if chrom not in chromosomes:
                    continue

            sorted_baited_digest_keys = self.get_baited_digest_keys_sorted_by_sta_pos()
            for baited_digest_key in sorted_baited_digest_keys[chrom]:
                baited_digest = self._baited_digest_dict[chrom][baited_digest_key]

                # Generate random numbers
                r_number_r = np.random.randint(155, 255)  # 100 colors
                r_number_b = np.random.randint(0, 100)

                # Directed interactions NE
                red_color = r_number_r
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_b
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('DI', 'NE'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['DI'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                # Directed interactions EN
                red_color = r_number_b
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_r
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('DI', 'EN'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['DI'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                # Undirected reference interactions NE
                red_color = r_number_r
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_b
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('UIR', 'NE'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['UIR'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                # Undirected reference interactions EN
                red_color = r_number_b
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_r
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('UIR', 'EN'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['UIR'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                # ALL interactions NE
                red_color = r_number_r
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_b
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('ALL', 'NE'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['ALL'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                # ALL interactions EN
                red_color = r_number_b
                green_color = 355 - r_number_r - r_number_b
                blue_color = r_number_r
                strand = "+"
                for d_inter in baited_digest.get_interactions_sorted_by_dist('ALL', 'EN'):
                    sta = str(d_inter.fromA)
                    end = str(d_inter.toB)
                    simple_twisted = str(d_inter.n_simple) + ':' + str(d_inter.n_twisted)
                    bed_file_streams['ALL'][bed_file_number].write(
                        chrom + '\t' + sta + '\t' + end + '\t' + simple_twisted + "\t0\t" + strand + "\t" + sta + '\t' + end +
                        '\t' + str(red_color) + "," + str(green_color) + "," + str(blue_color) + '\n')

                bed_file_number += 1
                if bed_file_number == 10:
                    bed_file_number = 1

        # Close streams for BED files
        for i in range(0, 10):
            bed_file_streams['DI'][i].close()
            bed_file_streams['UIR'][i].close()
            bed_file_streams['UI'][i].close()
            bed_file_streams['ALL'][i].close()

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
            proportion_di = round(n_di_interactions / float(n_total_interactions), 2)
            proportion_uir = round(n_uir_interactions / float(n_total_interactions), 2)
            di_ne_proportion = round(n_di_ne_interactions / float(n_di_interactions), 2)
            uir_ne_proportion = round(n_uir_ne_interactions / float(n_uir_interactions), 2)
            print(chr + '\t' +
                  str(n_total_interactions) + '\t' +
                  str(n_di_interactions) + '\t' +
                  str(proportion_di) + '\t' +
                  str(proportion_uir) + '\t' +
                  str(di_ne_proportion) + '\t' +
                  str(uir_ne_proportion)
                  )
