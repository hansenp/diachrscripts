from .diachromatic_interaction_set import DiachromaticInteractionSet
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import gzip


class IaFreqDistAnalysis:

    def __init__(self):

        # Define interaction and enrichment categories
        self.i_cats = [
            'DIX',
            'DI',
            'UIR',
            'UI',
            'ALL']
        self.i_cat_colors = [
            'red',
            'orange',
            'lightblue',
            'lightgray',
            'cornflowerblue']
        self.i_cat_names = [
            'Unbalanced without reference',
            'Unbalanced',
            'Balanced reference',
            'Balanced',
            'All']
        self.e_cats = [
            'NN',
            'NE',
            'EN',
            'EE']

        # Dictionary that contains ingested interaction grouped by chromosomes
        self._grouped_interactions = defaultdict()

        # Dictionary with information about ingestion of interactions
        self._ingest_interaction_set_info_dict = {
            'TOTAL_INTERACTIONS_READ': 0,
            'DIX': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
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
            }
        }

        # Prepare dictionary with read pair numbers used to create plots
        self.rp_num_dict = dict()

        # Prepare dictionary with interaction distances used to create plots
        self.i_dist_dict = dict()

    def ingest_interaction_set(self, d11_inter_set: DiachromaticInteractionSet, verbose: bool = False):
        """
        Ingests interactions from a 'DiachromaticInteractionSet' and groups them by chromosomes as well as interaction
        and enrichment category.

        :param d11_inter_set: DiachromaticInteractionSet with DiachromaticInteraction11 interactions
        :param verbose: If true, progress information will be displayed
        :return: Dictionary containing information about ingested  interactions
        """

        if verbose:
            print(
                "[INFO] Reading interactions and group them according to chromosomes, interaction and enrichment "
                "category ...")

        for d11_inter in d11_inter_set.interaction_list:
            self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] += 1
            if verbose:
                if self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] % 1000000 == 0:
                    print("\t[INFO] Read " + "{:,}".format(
                        self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + " interactions ...")

            # Get chromosome as well as interaction and enrichment category
            chrom = d11_inter.chrA
            i_cat_type = d11_inter.get_category()
            e_cat_type = d11_inter.enrichment_status_tag_pair

            # Count interaction type
            self._ingest_interaction_set_info_dict[i_cat_type][e_cat_type] += 1
            self._ingest_interaction_set_info_dict['ALL'][e_cat_type] += 1

            # Create a new dictionary for this chromosome, if this is the first interaction seen on this chromosome
            if d11_inter.chrA not in self._grouped_interactions:
                self._grouped_interactions[chrom] = defaultdict()
                for i_cat in self.i_cats:
                    self._grouped_interactions[chrom][i_cat] = dict()
                    for e_cat in self.e_cats:
                        self._grouped_interactions[chrom][i_cat][e_cat] = []

            # Add interaction to grouped interactions
            self._grouped_interactions[chrom][i_cat_type][e_cat_type].append(d11_inter)
            self._grouped_interactions[chrom]['ALL'][e_cat_type].append(d11_inter)

        if verbose:
            print("\t[INFO] Total number of interactions read: " + "{:,}".format(
                self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']))
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
        report += "\t[INFO] Broken down by interaction category and enrichment status: " + '\n'
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', ]:
            report += "\t\t[INFO] " + i_cat + ": " + '\n'
            for e_cat in ['NN', 'EE', 'NE', 'EN']:
                report += "\t\t\t[INFO] " + e_cat + ": " + "{:,}".format(self._ingest_interaction_set_info_dict[i_cat][e_cat]) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def _get_empty_num_dict(self,
                            i_cats: list = None,
                            i_cat_colors: list = None,
                            i_cat_names: list = None,
                            e_cats=None):
        """
        This function initializes a data structure that can be filled with interaction numbers or read pair numbers.

        The data structure consists of a dictionary with further sub-dictionaries. At the top level there are six
        dictionaries, one for each interaction category. In the level below, there are four lists for each category,
        one for each enrichment category (NN, NE, EN, EE).

        In addition to the lists with the pairs, there is a list with chromosomes that were taken into account.

        :return: A dictionary containing lists of numbers for the various interaction and enrichment categories and
        a list of chromosomes that were taken into account.
        """

        num_dict = dict()

        for i in range(len(i_cats)):
            num_dict[i_cats[i]] = dict()
            for e_cat in e_cats:
                num_dict[i_cats[i]][e_cat] = []
            if i_cat_names is not None:
                num_dict[i_cats[i]]['NAME'] = i_cat_names[i]
            if i_cat_colors is not None:
                num_dict[i_cats[i]]['COLOR'] = i_cat_colors[i]
        num_dict['CHROMOSOMES'] = []

        return num_dict

    def get_all_rp_nums_and_i_dists(self,
                                    chromosomes: [str] = None,
                                    verbose: bool = False):
        """
        This function collects read pair numbers and interaction distances in the various interaction and enrichment
        categories.

        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21']
        :param verbose: If true, messages about progress will be written to the screen
        :return: A dictionary containing the results and a list of chromosomes that were taken into account.
        """

        if verbose:
            print("[INFO] Getting all read pair numbers and interaction distances ...")

        # Reset num_dicts
        self.rp_num_dict = self._get_empty_num_dict(
            i_cats=self.i_cats,
            i_cat_colors=self.i_cat_colors,
            i_cat_names=self.i_cat_names,
            e_cats=self.e_cats)
        self.rp_num_dict['NUM_TYPE'] = 'Read pair number'

        # Reset num_dicts
        self.i_dist_dict = self._get_empty_num_dict(
            i_cats=self.i_cats,
            i_cat_colors=self.i_cat_colors,
            i_cat_names=self.i_cat_names,
            e_cats=self.e_cats)
        self.i_dist_dict['NUM_TYPE'] = 'Interaction distance'

        # Combine lists of read pair numbers and distances from all chromosomes
        for chrom in self._grouped_interactions.keys():
            if chromosomes is not None:
                if chrom not in chromosomes:
                    continue

            if verbose:
                print("\t[INFO] Processing chromosome " + chrom + " ...")

            self.rp_num_dict['CHROMOSOMES'].append(chrom)
            self.i_dist_dict['CHROMOSOMES'].append(chrom)
            for i_cat in ['DIX', 'DI','UIR', 'UI', 'ALL']:
                for e_cat in ['NN', 'NE', 'EN', 'EE']:
                    for d11_inter in self._grouped_interactions[chrom][i_cat][e_cat]:
                        self.rp_num_dict[i_cat][e_cat].append(d11_inter.rp_total)
                        self.i_dist_dict[i_cat][e_cat].append(d11_inter.i_dist)

        if verbose:
            print("[INFO] ... done.")

        return self.rp_num_dict, self.i_dist_dict

    def get_all_rp_nums_or_i_dists_histograms(self,
                                              num_dict: dict = None,
                                              i_cats: list = None,
                                              e_cats: list = None,
                                              q_lim: float = 0.95,
                                              description: str = "DESCRIPTION",
                                              pdf_file_name: str = "num_histograms.pdf"):
        """
        This function creates the histograms for read pair numbers or interaction distances in the various
        interaction and enrichment categories.

        :param i_cats: List of interaction categories
        :param e_cats: List of enrichment categories
        :param num_dict: Dictionary that was created with the function 'get_all_rp_nums_or_i_dists()'
        :param q_lim: Chose upper limit for x-axes based on quantile (affects presentation only)
        :param description: Brief description that is shown in the plot above the histograms
        :param pdf_file_name: Name of the PDF file that will be created.
        :return: A matplotlib 'Figure' object that can be displayed in Jupyter notebooks
        """

        print("1")
        # Catch wrong input
        allowed_i_cats = ['DIX', 'DI', 'UI', 'UIR', 'ALL']
        for i_cat in i_cats:
            if i_cat not in allowed_i_cats:
                print("[ERROR] Illegal interaction category tag! Allowed: 'DIX', 'DI', 'D_S', 'D_T', 'UI', 'UIR' and 'ALL'")
                return

        allowed_e_cats = ['NN', 'NE', 'EN', 'EE']
        for e_cat in e_cats:
            if e_cat not in allowed_e_cats:
                print("[ERROR] Illegal interaction enrichment tag! Allowed: 'NN', 'NE', 'EN', 'EE'")
                return

        print("2")
        # Prepare grid for individual plots
        n = len(i_cats)
        y = 11.79/4.75
        fig_height = 0.75 * y + n * y
        m = len(e_cats)
        x = 4.75
        fig_width = m * x
        fig, ax = plt.subplots(nrows=n+1, ncols=m, figsize=(fig_width, fig_height),
                               gridspec_kw={'height_ratios': [0.75] + [1]*n})

        # Hide unnecessary subplots
        for j in range(0, m):
            ax[0][j].axis('off')

        # Determine bin size
        x_lim = 0
        for i in range(0, n):
            for j in range(0, m):
                q = np.quantile(num_dict[i_cats[i]][e_cats[j]], q_lim)
                if x_lim < q:
                    x_lim = q
        x_ticks, x_tick_labels = self.make_ticks(x_lim)
        bin_width = int(x_lim / 30)
        if bin_width < 1:
            bin_width = 1

        print("3")

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

        fig.text(0.015, 1-((1-0.97)*(11.79/fig_height)), num_dict['NUM_TYPE'] + 's', fontsize=18,
                 fontweight='bold')
        fig.text(0.030, 1-((1-0.94)*(11.79/fig_height)), 'Description: ' + description, fontsize=12)
        fig.text(0.030, 1-((1-0.92)*(11.79/fig_height)), 'Bin size: ' + "{:,}".format(bin_width), fontsize=12)
        fig.text(0.030, 1-((1-0.90)*(11.79/fig_height)), 'For chromosomes:', fontsize=12)
        if len(num_dict['CHROMOSOMES']) < 22:
            fig.text(0.045, 1-((1-0.88)*(11.79/fig_height)), '[' + ", ".join(i for i in num_dict['CHROMOSOMES']) + ']', fontsize=8)
        else:
            fig.text(0.045, 1-((1-0.88)*(11.79/fig_height)), '[' + ", ".join(i for i in num_dict['CHROMOSOMES'][:22]) + ',', fontsize=8)
            fig.text(0.045, 1-((1-0.86)*(11.79/fig_height)), ", ".join(i for i in num_dict['CHROMOSOMES'][22:]) + ']', fontsize=8)

        # Create histograms for each interaction category
        # -----------------------------------------------
        print("4")
        # Prepare bins
        x_max = 0
        for i in range(0, n):
            for j in range(0, m):
                if x_max < max(num_dict[i_cats[i]][e_cats[j]]):
                    x_max = max(num_dict[i_cats[i]][e_cats[j]])
        bins = range(0, x_max + bin_width, bin_width)

        abs_2d_array = [[] for i in range(0, n)]
        densities_2d_array = [[] for i in range(0, n)]
        for i in range(0, n):

            print("4, i=" + str(i))

            for j in range(0, m):
                print("4, j=" + str(j))

                # Create histogram
                counts, bins, patches = ax[i + 1][j].hist(
                    num_dict[i_cats[i]][e_cats[j]],
                    bins=bins, density=False,
                    facecolor=num_dict[i_cats[i]]['COLOR'],
                    edgecolor="dimgray",
                    linewidth=0.5,
                    alpha=1)
                ax[i + 1][j].set_title(
                    num_dict[i_cats[i]]['NAME'] + ' - ' + e_cats[j],
                    loc='left')
                ax[i + 1][j].set_xlabel(num_dict['NUM_TYPE'])
                ax[i + 1][j].set_ylabel('Frequency')
                ax[i + 1][j].set_xticks(x_ticks)
                ax[i + 1][j].set_xticklabels(x_tick_labels)
                ax[i + 1][j].set_xlim(0, x_lim)

                # Draw vertical lines and shaded areas for median and MAD
                median = np.median(num_dict[i_cats[i]][e_cats[j]])
                mad = stats.median_abs_deviation(num_dict[i_cats[i]][e_cats[j]])
                ax[i + 1][j].axvline(median, linestyle='--', linewidth=0.75, color='blue', zorder=2)
                ax[i + 1][j].axvspan(median, median + mad, color='green', alpha=0.25, zorder=0)

                # Keep track of bin counts and densities
                abs_2d_array[i].append(counts)
                sum_counts = sum(counts)
                densities_2d_array[i].append([count / sum_counts for count in counts]) # slow?

        print("5")
        # Add second axes with densities and normalize all histograms to maximum density
        # ------------------------------------------------------------------------------

        # Transform bins to center bin positions. The original list for the bins contains the boundaries between bins,
        # but for the bar plot we need the center positions of bins.
        bcp = bin_width / 2
        bcp_list = [bcp]
        for i in range(len(densities_2d_array[0][0]) - 1):
            bcp += bin_width
            bcp_list.append(bcp)

        # Determine maximal density in all histograms
        yd_max = 0
        for i in range(0, n):
            for j in range(0, m):
                if yd_max < max(densities_2d_array[i][j]):
                    yd_max = max(densities_2d_array[i][j])
        y_padding = yd_max / 20
        yd_max += y_padding

        print("6")

        # Add density axes, normalize and add text labels to histograms
        for i in range(0, n):
            for j in range(0, m):

                yc_max = yd_max * sum(abs_2d_array[i][j])
                yc_ticks, yc_tick_labels = self.make_ticks(yc_max)
                ax[i + 1][j].set_yticks(yc_ticks)
                ax[i + 1][j].set_yticklabels(yc_tick_labels)
                ax[i + 1][j].set_ylim(0, yc_max)
                ax_dens = ax[i + 1][j].twinx()
                #ax_dens.plot(bcp_list, densities_2d_array[i][j], color='red', linewidth=0.5)
                ax_dens.set_ylabel('Density', labelpad=7)
                ax_dens.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
                ax_dens.set_ylim(0, yd_max)
                # Add text labels with total read pair or interaction numbers, median and median absolute deviation
                ax[i + 1][j].text(x_lim - (x_lim / 3),
                                  yc_max - (yc_max / 3.2), #3.2
                                  'n: ' + "{:,}".format(len(num_dict[i_cats[i]][e_cats[j]])) + '\n' +
                                  'Mdn: ' + "{:,.0f}".format(np.median(num_dict[i_cats[i]][e_cats[j]])) + '\n' +
                                  'MAD: ' + "{:,.0f}".format(stats.median_abs_deviation(num_dict[i_cats[i]][e_cats[j]])),
                                  fontsize=8.5,
                                  bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        print("7")
        # Save and return figure
        fig.tight_layout(pad=1.25)
        fig.savefig(pdf_file_name)
        return fig

    def get_all_rp_nums_or_i_dists_denisty_diff_plot(self,
                                                     num_dict: dict = None,
                                                     i_cats: list = None,
                                                     e_cats: list = None,
                                                     q_lim: float = 0.95,
                                                     description: str = "DESCRIPTION",
                                                     pdf_file_name: str = "density_diff_plot.pdf"):
        """
        This function creates the density difference plot for given pairs of interaction and enrichment categories.

        :param num_dict: A data structure that contains read pair numbers and distances for all interaction categories
        :param i_cats: List of two interaction categories to be compared, e.g. ['DI', 'UIR']
        :param e_cats: List of two enrichment categories to be compared, e.g. ['EE', 'EN']
        :param q_lim: Upper limit at quantile level, e.g. at 0.95, 5% of the data with the larger quantile are not shown
        :param description: Brief description that is shown in the plot above the histograms
        :param pdf_file_name: Name of the PDF file that will be created
        :return: A 'Figure' object of 'matplotlib' that can be displayed in a Jupyter notebook
        """

        # Catch wrong input
        if len(i_cats) !=2:
            print("[ERROR] The density difference plot is only defined for two interaction categories!")
            return
        if len(e_cats) !=2:
            print("[ERROR] The density difference plot is only defined for two enrichment categories!")
            return

        print("1")

        # Prepare grid for individual plots
        fig_height = 9.16
        n = len(i_cats)
        m = len(e_cats)
        fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(14.25, fig_height),
                               gridspec_kw={'height_ratios': [0.75, 1, 1, 1]})

        # Hide unnecessary subplots
        ax[0][2].axis('off')
        ax[3][2].axis('off')

        # Determine bin size
        x_lim = 0
        for i in range(0, n):
            for j in range(0, m):
                q = np.quantile(num_dict[i_cats[i]][e_cats[j]], q_lim)
                if x_lim < q:
                    x_lim = q
        x_ticks, x_tick_labels = self.make_ticks(x_lim)
        bin_width = int(x_lim / 30)

        print("2")

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
        fig.text(0.015*(9.5/14.25), 1-((1-0.97)*(11.79/fig_height)), num_dict['NUM_TYPE'] + 's', fontsize=18,
                 fontweight='bold')
        fig.text(0.030*(9.5/14.25), 1-((1-0.94)*(11.79/fig_height)), 'Description: ' + description, fontsize=12)
        fig.text(0.030 * (9.5 / 14.25), 1-((1-0.92)*(11.79/fig_height)), 'Bin size: ' + "{:,}".format(bin_width), fontsize=12)
        fig.text(0.030*(9.5/14.25), 1-((1-0.90)*(11.79/fig_height)), 'For chromosomes:', fontsize=12)
        fig.text(0.045*(9.5/14.25), 1-((1-0.88)*(11.79/fig_height)), '[' + ", ".join(i for i in num_dict['CHROMOSOMES']) + ']', fontsize=8)

        print("3")

        # Prepare bins
        x_max = 0
        for i in range(0, n):
            for j in range(0, m):
                if x_max < max(num_dict[i_cats[i]][e_cats[j]]):
                    x_max = max(num_dict[i_cats[i]][e_cats[j]])
        bins = range(0, x_max + bin_width, bin_width)


        # Create histograms for the two categories
        # ----------------------------------------

        print("4")

        abs_2d_array = [[] for i in range(0, n)]
        densities_2d_array = [[] for i in range(0, n)]
        for i in range(0, n):
            print("4" + ': i=' + str(i))

            for j in range(0, m):
                print("4" + ': i=' + str(i) + ', j=' + str(j))

                # Create histogram
                counts, bins, patches = ax[i + 1][j].hist(
                    num_dict[i_cats[i]][e_cats[j]],
                    bins=bins, density=False,
                    facecolor=num_dict[i_cats[i]]['COLOR'],
                    edgecolor="dimgray",
                    linewidth=0.5,
                    alpha=1)
                ax[i + 1][j].set_title(
                    num_dict[i_cats[i]]['NAME'] + ' - ' + e_cats[j],
                    loc='left')
                ax[i + 1][j].set_xlabel(num_dict['NUM_TYPE'])
                ax[i + 1][j].set_ylabel('Frequency')
                ax[i + 1][j].set_xticks(x_ticks)
                ax[i + 1][j].set_xticklabels(x_tick_labels)
                ax[i + 1][j].set_xlim(0, x_lim)

                # Draw vertical lines and shaded areas for median and MAD
                median = np.median(num_dict[i_cats[i]][e_cats[j]])
                mad = stats.median_abs_deviation(num_dict[i_cats[i]][e_cats[j]])
                ax[i + 1][j].axvline(median, linestyle='--', linewidth=0.75, color='blue', zorder=2)
                ax[i + 1][j].axvspan(median, median + mad, color='green', alpha=0.25, zorder=0)

                # Keep track of bin counts and densities
                abs_2d_array[i].append(counts)
                sum_counts = sum(counts)
                densities_2d_array[i].append([count / sum_counts for count in counts])

        print("5")

        # Add second axes with densities and normalize all histograms to maximum density
        # ------------------------------------------------------------------------------

        # Transform bins to center bin positions. The original list for the bins contains the boundaries between bins,
        # but for the bar plot we need the center positions of bins.
        bcp = bin_width / 2
        bcp_list = [bcp]
        for i in range(len(densities_2d_array[0][0]) - 1):
            bcp += bin_width
            bcp_list.append(bcp)

        print("6")

        # Determine maximal density in all four plots
        yd_max = 0
        for i in range(0, n):
            for j in range(0, m):
                if yd_max < max(densities_2d_array[i][j]):
                    yd_max = max(densities_2d_array[i][j])
        y_padding = yd_max / 20
        yd_max += y_padding

        print("7")

        # Add density axes, normalize and add text labels to histograms
        for i in range(0, n):
            for j in range(0, m):
                yc_max = yd_max * sum(abs_2d_array[i][j])
                yc_ticks, yc_tick_labels = self.make_ticks(yc_max)
                ax[i + 1][j].set_yticks(yc_ticks)
                ax[i + 1][j].set_yticklabels(yc_tick_labels)
                ax[i + 1][j].set_ylim(0, yc_max)
                ax_dens = ax[i + 1][j].twinx()
                #ax_dens.plot(bcp_list, densities_2d_array[i][j], color='red', linewidth=0.5)
                ax_dens.set_ylabel('Density', labelpad=7)
                ax_dens.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
                ax_dens.set_ylim(0, yd_max)
                # Add text labels with total read pair or interaction numbers, median and median absolute deviation
                ax[i + 1][j].text(x_lim - (x_lim / 3),
                                  yc_max - (yc_max / 3.2), #3.2
                                  'n: ' + "{:,}".format(len(num_dict[i_cats[i]][e_cats[j]])) + '\n' +
                                  #'Q1: ' + "{:,.0f}".format(np.quantile(num_dict[i_cats[i]][e_cats[j]], 0.25)) + '\n' +
                                  'Mdn: ' + "{:,.0f}".format(np.median(num_dict[i_cats[i]][e_cats[j]])) + '\n' +
                                  'MAD: ' + "{:,.0f}".format(stats.median_abs_deviation(num_dict[i_cats[i]][e_cats[j]])),
                                  fontsize=8.5,
                                  bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        # Add subplots for density differences
        # ------------------------------------
        print("8")

        # Get four lists with densities differences from all four density lists
        density_diff_cat_0_cat_1_ne = [] # Below histograms on the left
        density_diff_cat_0_cat_1_en = [] # Below histograms on the right
        density_diff_cat_0_ne_en = []    # To the left of the upper histograms
        density_diff_cat_1_ne_en = []    # To the left of the lower histograms
        for i in range(0, len(densities_2d_array[0][0])):
            density_diff_cat_0_cat_1_ne.append(densities_2d_array[0][0][i] - densities_2d_array[1][0][i])
            density_diff_cat_0_cat_1_en.append(densities_2d_array[0][1][i] - densities_2d_array[1][1][i])
            density_diff_cat_0_ne_en.append(densities_2d_array[0][0][i] - densities_2d_array[0][1][i])
            density_diff_cat_1_ne_en.append(densities_2d_array[1][0][i] - densities_2d_array[1][1][i])

        print("9")

        # Perform KS test
        print(i_cats)
        print(e_cats)
        print(stats.ks_2samp(num_dict[i_cats[0]][e_cats[0]], num_dict[i_cats[1]][e_cats[0]]))
        print(stats.ks_2samp(num_dict[i_cats[0]][e_cats[1]], num_dict[i_cats[1]][e_cats[1]]))
        print(stats.ks_2samp(num_dict[i_cats[0]][e_cats[0]], num_dict[i_cats[0]][e_cats[1]]))
        print(stats.ks_2samp(num_dict[i_cats[1]][e_cats[0]], num_dict[i_cats[1]][e_cats[1]]))

        print("10")

        # Determine sum of density differences
        dd_cat_0_cat_1_ne_sum = sum(map(abs, density_diff_cat_0_cat_1_ne))
        dd_cat_0_cat_1_en_sum = sum(map(abs, density_diff_cat_0_cat_1_en))
        dd_cat_0_ne_en_sum = sum(map(abs, density_diff_cat_0_ne_en))
        dd_cat_1_ne_en_sum = sum(map(abs, density_diff_cat_1_ne_en))

        print("11")

        # Bar plot below histograms on the left
        ax[3][0].bar(bcp_list,
                     density_diff_cat_0_cat_1_ne,
                     width=bin_width,
                     color=[num_dict[i_cats[0]]['COLOR'] if 0 < dd else num_dict[i_cats[1]]['COLOR'] for dd in density_diff_cat_0_cat_1_ne],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[3][0].set_title(e_cats[0], loc='left')
        ax[3][0].set_xticks(x_ticks)
        ax[3][0].set_xticklabels(x_tick_labels)
        ax[3][0].set_xlim(0, x_lim)
        ax[3][0].set_xlabel(num_dict['NUM_TYPE'])
        ax[3][0].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[3][0].ticklabel_format(axis='y', style='sci', scilimits=(-2,0))
        ax[3][0].yaxis.tick_right()
        ax[3][0].yaxis.set_ticks_position('both')
        ax[3][0].set_ylabel('Density difference', labelpad=7)
        ax[3][0].yaxis.set_label_position('right')

        # Bar plot below histograms on the right
        ax[3][1].bar(bcp_list,
                     density_diff_cat_0_cat_1_en,
                     width=bin_width,
                     color=[num_dict[i_cats[0]]['COLOR'] if 0 < dd else num_dict[i_cats[1]]['COLOR'] for dd in density_diff_cat_0_cat_1_en],
                     edgecolor="dimgray",
                     linewidth=0.5
                     )
        ax[3][1].set_title(e_cats[1], loc='left')
        ax[3][1].set_xticks(x_ticks)
        ax[3][1].set_xticklabels(x_tick_labels)
        ax[3][1].set_xlim(0, x_lim)
        ax[3][1].set_xlabel(num_dict['NUM_TYPE'])
        ax[3][1].axhline(0, linestyle='-.', linewidth=0.75, color='blue', zorder=2)
        ax[3][1].ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
        ax[3][1].yaxis.tick_right()
        ax[3][1].yaxis.set_ticks_position('both')
        ax[3][1].set_ylabel('Density difference', labelpad=7)
        ax[3][1].yaxis.set_label_position('right')

        print("12")

        # Make y-axes comparable for the two bar plots below the histograms and add sums of density differences
        y_min = min(min(density_diff_cat_0_cat_1_ne), min(density_diff_cat_0_cat_1_en))
        y_max = max(max(density_diff_cat_0_cat_1_ne), max(density_diff_cat_0_cat_1_en))
        y_padding = ((y_max - y_min)/20)
        y_min -= y_padding
        y_max += y_padding
        ax[3][0].set_ylim(y_min, y_max)
        ax[3][1].set_ylim(y_min, y_max)
        #ax[3][0].text(x_lim - (x_lim / 3),
        #              y_max - ((y_max - y_min) / 8),
        #              'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_cat_1_ne_sum),
        #              fontsize=9,
        #              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
        #ax[3][1].text(x_lim - (x_lim / 3),
        #              y_max - ((y_max - y_min) / 8),
        #              'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_cat_1_en_sum),
        #              fontsize=9,
        #              bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

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
        ax[1][2].ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
        ax[1][2].yaxis.tick_right()
        ax[1][2].yaxis.set_ticks_position('both')
        ax[1][2].set_ylabel('Density difference', labelpad=7)
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
        ax[2][2].ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
        ax[2][2].yaxis.tick_right()
        ax[2][2].yaxis.set_ticks_position('both')
        ax[2][2].set_ylabel('Density difference', labelpad=7)
        ax[2][2].yaxis.set_label_position('right')

        print("13")

        # Make y-axes comparable for the two bar plots to the left of the histograms and add sums of density differences
        y_min = min(min(density_diff_cat_0_ne_en), min(density_diff_cat_1_ne_en))
        y_max = max(max(density_diff_cat_0_ne_en), max(density_diff_cat_1_ne_en))
        y_padding = ((y_max - y_min) / 20)
        y_min -= y_padding
        y_max += y_padding
        ax[1][2].set_ylim(y_min, y_max)
        ax[2][2].set_ylim(y_min, y_max)
        # ax[1][2].text(x_lim - (x_lim / 3),
        #               y_max - ((y_max - y_min) / 8),
        #               'sum(|dd|): ' + "{:.2f}".format(dd_cat_0_ne_en_sum),
        #               fontsize=9,
        #               bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))
        # ax[2][2].text(x_lim - (x_lim / 3),
        #               y_max - ((y_max - y_min) / 8),
        #               'sum(|dd|): ' + "{:.2f}".format(dd_cat_1_ne_en_sum),
        #               fontsize=9,
        #               bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, boxstyle='round'))

        print("14")

        # Save and return figure
        fig.tight_layout(pad=1.25)
        fig.savefig(pdf_file_name)
        return fig

    def write_num_table(self,
                        out_prefix: str = 'OUT_PREFIX',
                        description: str =  'DESCRIPTION',
                        chromosomes: [str] = None,
                        verbose: bool = False):
        """

        :param out_prefix: Two tab separated files with this prefix will be created, one for read pair numbers and
        one for interaction distances
        :param description: A short description tag, e.g. 'MK'
        :param chromosomes: The analysis can be restricted to subsets chromosomes, e.g. ['chr19', 'chr20', 'chr21']
        :param verbose: If true, messages about progress will be written to the screen
        :return: Nothing
        """

        if verbose:
            print("[INFO] Writing results about read pair numbers and interaction distances to tab separated files ...")

        #  Create header row
        header_row = ''
        for i_cat in self.i_cats:
            for e_cat in self.e_cats:
                header_row += '\t' + i_cat + '|' + e_cat + '|N'
                header_row += '\t' + i_cat + '|' + e_cat + '|MED'
                header_row += '\t' + i_cat + '|' + e_cat + '|MAD'
        header_row += '\t' + 'CHROMS'

        # Create two files, one for read pair numbers and one for interaction distances
        out_rp_num_stats = open(out_prefix + "_rp_num_stats.txt", 'wt')
        out_rp_num_stats.write('DESCRIPTION' + '\t' + 'DESCRIPTION_SHORT' + header_row  + '\n')
        out_i_dist_stats = open(out_prefix + "_i_dist_stats.txt", 'wt')
        out_i_dist_stats.write('DESCRIPTION' + '\t' + 'DESCRIPTION_SHORT' + header_row  + '\n')

        # Create one read pair number num_dict for all chromosomes
        rp_num_dict = self._get_empty_num_dict(
            i_cats=self.i_cats,
            i_cat_colors=self.i_cat_colors,
            i_cat_names=self.i_cat_names,
            e_cats=self.e_cats)

        # Create one read pair number num_dict for all chromosomes
        i_dist_dict = self._get_empty_num_dict(
            i_cats=self.i_cats,
            i_cat_colors=self.i_cat_colors,
            i_cat_names=self.i_cat_names,
            e_cats=self.e_cats)

        for chrom in self._grouped_interactions.keys():
            if chromosomes is not None:
                if chrom not in chromosomes:
                    continue

            if verbose:
                print("\t[INFO] Processing chromosome " + chrom + " ...")

            # Create one read pair number num_dict for all chromosomes
            rp_num_dict_chr = self._get_empty_num_dict(
                i_cats=self.i_cats,
                i_cat_colors=self.i_cat_colors,
                i_cat_names=self.i_cats,
                e_cats=self.e_cats)
            rp_num_dict_chr['NUM_TYPE'] = 'Read pair number'
            rp_num_dict_chr['CHROMOSOMES'].append(chrom)

            # Create one read pair number num_dict for all chromosomes
            i_dist_dict_chr = self._get_empty_num_dict(
                i_cats=self.i_cats,
                i_cat_colors=self.i_cat_colors,
                i_cat_names=self.i_cats,
                e_cats=self.e_cats)
            i_dist_dict_chr['NUM_TYPE'] = 'Interaction distance'
            i_dist_dict_chr['CHROMOSOMES'].append(chrom)

            rp_num_dict['CHROMOSOMES'].append(chrom)
            i_dist_dict['CHROMOSOMES'].append(chrom)
            for i_cat in self.i_cats:
                for e_cat in self.e_cats:
                    for d11_inter in self._grouped_interactions[chrom][i_cat][e_cat]:
                        rp_num_dict_chr[i_cat][e_cat].append(d11_inter.rp_total)
                        i_dist_dict_chr[i_cat][e_cat].append(d11_inter.i_dist)
                        rp_num_dict[i_cat][e_cat].append(d11_inter.rp_total)
                        i_dist_dict[i_cat][e_cat].append(d11_inter.i_dist)

            # Get table row for num_dicts of this chromosome
            rp_num_chr_table_row = self.get_num_dict_table_row(description = description + '|RP_NUM', description_short = description, num_dict = rp_num_dict_chr)
            out_rp_num_stats.write(rp_num_chr_table_row + '\n')
            i_dist_num_chr_table_row = self.get_num_dict_table_row(description = description + '|I_DIST', description_short = description, num_dict = i_dist_dict_chr)
            out_i_dist_stats.write(i_dist_num_chr_table_row + '\n')

        # Get table row for num_dict of all chromosomes
        rp_num_table_row = self.get_num_dict_table_row(description = description + '|RP_NUM', description_short = description, num_dict = rp_num_dict)
        out_rp_num_stats.write(rp_num_table_row + '\n')
        i_dist_num_table_row = self.get_num_dict_table_row(description = description + '|I_DIST', description_short = description, num_dict = i_dist_dict)
        out_i_dist_stats.write(i_dist_num_table_row + '\n')

        # Close the two files with tables
        out_rp_num_stats.close()
        out_i_dist_stats.close()

        if verbose:
            print("\t[INFO] Generated files:")
            print("\t[INFO] For read pair numbers:")
            print("\t   " + out_prefix + "_rp_num_stats.txt")
            print("\t[INFO] For interaction distances:")
            print("\t   " + out_prefix + "_i_dist_stats.txt")
            print("[INFO] ... done.")

    def get_num_dict_table_row(self,
                               description: str =  'DESCRIPTION',
                               description_short: str =  'CT',
                               num_dict = dict):

        # Create description tag for first column
        if len(num_dict['CHROMOSOMES']) == 1:
            table_row = description + '|' + num_dict['CHROMOSOMES'][0].upper() + '\t'
        else:
            table_row = description + '|CHR_ALL'  + '\t'

        table_row += description_short  + '\t'

        # Add values for all interaction and enrichment categories
        for i_cat in self.i_cats:
            for e_cat in self.e_cats:
                n = str(len(num_dict[i_cat][e_cat]))
                if 0 < int(n):
                    median = "{:.0f}".format(np.median(num_dict[i_cat][e_cat]))
                    mad = "{:.0f}".format(stats.median_abs_deviation(num_dict[i_cat][e_cat]))
                else:
                    median = 'NA'
                    mad = 'NA'
                table_row += n + '\t' +  median + '\t' +  mad + '\t'

        # Add last column for chromosomes that were taken into aaccount
        table_row += str(num_dict['CHROMOSOMES']).replace(' ', '')

        return table_row

    def read_num_table(self):
        pass

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
