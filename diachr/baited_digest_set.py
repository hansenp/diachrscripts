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

        if (number_pair_type is None) or number_pair_type not in ['I_NUM', 'RP_NUM', 'MED_RP_NUM', 'MED_I_DIST', 'DIG_LEN']:
            print("[ERROR] Invalid number pair type! Use one of the following:")
            print("\t[ERROR] 'I_NUM' - Interaction numbers")
            print("\t[ERROR] 'RP_NUM' - Read pair numbers")
            print("\t[ERROR] 'MED_RP_NUM' - Median read pair numbers")
            print("\t[ERROR] 'MED_I_DIST' - Median interaction distances")
            print("\t[ERROR] 'DIG_LEN' - Digest length")
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
        if number_pair_type == 'DIG_LEN':
            pair_dict['NUM_PAIR_TYPE'] = 'Digest length'

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
                    if number_pair_type == 'DIG_LEN':
                        pair_dict[i_cat]['NE'].extend(baited_digest.get_dig_len_list(i_cat, 'D1'))
                        pair_dict[i_cat]['EN'].extend(baited_digest.get_dig_len_list(i_cat, 'D2'))

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

    def get_neen_pvals(self, i_cat):
        for chrom in self._baited_digest_dict.keys():
            for key, baited_digest in self._baited_digest_dict[chrom].items():
                ne_num, en_num, p_val = baited_digest.get_neen_p_value(i_cat)
                if(4 < p_val):
                    print(key + '\t' + str(ne_num) + '\t' + str(en_num) + '\t' + str(p_val))