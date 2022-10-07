import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from .diachromatic_interaction_set import DiachromaticInteractionSet


class ReadTypeAndConfigCounter:

    def __init__(self):
        self.x = 1

    def test(self):
        print('Hurz')

    @staticmethod
    def count_read_types(d11_interaction_set: DiachromaticInteractionSet = None):

        # Initialize dictionary for counts
        rp_type_freq_dict = dict()
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B']:
            rp_type_freq_dict[i_cat] = dict()
            for e_cat in ['NN', 'EE', 'NE', 'EN', 'ALL']:
                rp_type_freq_dict[i_cat][e_cat] = dict()
                for rp_type in ['T0', 'T1', 'T2', 'T3']:
                    rp_type_freq_dict[i_cat][e_cat][rp_type] = 0

        # Iterate over all interactions and count types
        for d11_inter in d11_interaction_set.interaction_list:

            # Get tags for interaction and enrichment category
            i_cat = d11_inter.get_category()
            e_cat = d11_inter.enrichment_status_tag_pair

            # Counts within each category
            rp_type_freq_dict[i_cat][e_cat]['T0'] += d11_inter._simple_1
            rp_type_freq_dict[i_cat][e_cat]['T1'] += d11_inter._simple_2
            rp_type_freq_dict[i_cat][e_cat]['T2'] += d11_inter._twisted_1
            rp_type_freq_dict[i_cat][e_cat]['T3'] += d11_inter._twisted_2

            # Counts combined for interaction categories
            rp_type_freq_dict['ALL'][e_cat]['T0'] += d11_inter._simple_1
            rp_type_freq_dict['ALL'][e_cat]['T1'] += d11_inter._simple_2
            rp_type_freq_dict['ALL'][e_cat]['T2'] += d11_inter._twisted_1
            rp_type_freq_dict['ALL'][e_cat]['T3'] += d11_inter._twisted_2

            # Counts combined for enrichment categories
            rp_type_freq_dict[i_cat]['ALL']['T0'] += d11_inter._simple_1
            rp_type_freq_dict[i_cat]['ALL']['T1'] += d11_inter._simple_2
            rp_type_freq_dict[i_cat]['ALL']['T2'] += d11_inter._twisted_1
            rp_type_freq_dict[i_cat]['ALL']['T3'] += d11_inter._twisted_2

            # Counts combined for interaction and enrichment categories
            rp_type_freq_dict['ALL']['ALL']['T0'] += d11_inter._simple_1
            rp_type_freq_dict['ALL']['ALL']['T1'] += d11_inter._simple_2
            rp_type_freq_dict['ALL']['ALL']['T2'] += d11_inter._twisted_1
            rp_type_freq_dict['ALL']['ALL']['T3'] += d11_inter._twisted_2

            # Counts for all unbalanced (DIX and DI)
            if i_cat == 'DIX' or i_cat == 'DI':
                rp_type_freq_dict['U'][e_cat]['T0'] += d11_inter._simple_1
                rp_type_freq_dict['U'][e_cat]['T1'] += d11_inter._simple_2
                rp_type_freq_dict['U'][e_cat]['T2'] += d11_inter._twisted_1
                rp_type_freq_dict['U'][e_cat]['T3'] += d11_inter._twisted_2
                rp_type_freq_dict['U']['ALL']['T0'] += d11_inter._simple_1
                rp_type_freq_dict['U']['ALL']['T1'] += d11_inter._simple_2
                rp_type_freq_dict['U']['ALL']['T2'] += d11_inter._twisted_1
                rp_type_freq_dict['U']['ALL']['T3'] += d11_inter._twisted_2

            # Counts for all balanced (UIR and UI) and balanced interactions
            if i_cat == 'UIR' or i_cat == 'UI':
                rp_type_freq_dict['B'][e_cat]['T0'] += d11_inter._simple_1
                rp_type_freq_dict['B'][e_cat]['T1'] += d11_inter._simple_2
                rp_type_freq_dict['B'][e_cat]['T2'] += d11_inter._twisted_1
                rp_type_freq_dict['B'][e_cat]['T3'] += d11_inter._twisted_2
                rp_type_freq_dict['B']['ALL']['T0'] += d11_inter._simple_1
                rp_type_freq_dict['B']['ALL']['T1'] += d11_inter._simple_2
                rp_type_freq_dict['B']['ALL']['T2'] += d11_inter._twisted_1
                rp_type_freq_dict['B']['ALL']['T3'] += d11_inter._twisted_2

        # Fill second dictionary with relative frequencies
        rp_type_dens_dict = copy.deepcopy(rp_type_freq_dict)
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B']:
            for e_cat in ['NN', 'EE', 'NE', 'EN', 'ALL']:
                rp_total = sum(rp_type_freq_dict[i_cat][e_cat].values())
                if 0 < rp_total:
                    rp_type_dens_dict[i_cat][e_cat]['T0'] = rp_type_freq_dict[i_cat][e_cat]['T0'] / rp_total
                    rp_type_dens_dict[i_cat][e_cat]['T1'] = rp_type_freq_dict[i_cat][e_cat]['T1'] / rp_total
                    rp_type_dens_dict[i_cat][e_cat]['T2'] = rp_type_freq_dict[i_cat][e_cat]['T2'] / rp_total
                    rp_type_dens_dict[i_cat][e_cat]['T3'] = rp_type_freq_dict[i_cat][e_cat]['T3'] / rp_total
                else:
                    rp_type_dens_dict[i_cat][e_cat]['T0'] = 0.0
                    rp_type_dens_dict[i_cat][e_cat]['T1'] = 0.0
                    rp_type_dens_dict[i_cat][e_cat]['T2'] = 0.0
                    rp_type_dens_dict[i_cat][e_cat]['T3'] = 0.0

        return rp_type_freq_dict, rp_type_dens_dict

    @staticmethod
    def print_read_type_frequency_tables(rp_type_freq_dict=None,
                                         rp_type_dens_dict=None,
                                         i_cats=['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B'],
                                         e_cats=['NN', 'EE', 'NE', 'EN', 'ALL'],
                                         print_dens=True):
        for i_cat in i_cats:
            print(i_cat)
            print('\t', end='')
            for e_cat in ['T0', 'T1', 'T2', 'T3']:
                print('\t' + e_cat, end='')
            print()
            for e_cat in e_cats:
                print('\t' + e_cat + ': ', end='')
                for ht_tag in ['T0', 'T1', 'T2', 'T3']:
                    if print_dens:
                        print('\t' + "{:.2f}".format(rp_type_dens_dict[i_cat][e_cat][ht_tag]), end='')
                    else:
                        print('\t' + "{:,}".format(rp_type_freq_dict[i_cat][e_cat][ht_tag]), end='')
                print()
            print()

    @staticmethod
    def create_read_type_frequency_bar_charts(
            rp_type_freq_dict=None,
            rp_type_dens_dict=None,
            i_cats=['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B'],
            e_cat='ALL',
            pdf_file_name='rp_type_frequency_bar_charts.pdf'):

        # Specify colors and labels
        rp_cat_colors = [(255 / 255, 160 / 255, 200 / 255), (255 / 255, 80 / 255, 120 / 255),
                         (80 / 255, 190 / 255, 120 / 255), (60 / 255, 150 / 255, 120 / 255)]
        rp_cat_labels = ['T0', 'T1', 'T2', 'T3']
        i_cat_titles = {'DIX': 'Unbalanced without reference',
                        'DI': 'Unbalanced with reference',
                        'UIR': 'Balanced reference',
                        'UI': 'Unbalanced no reference',
                        'ALL': 'All',
                        'U': 'Unbalanced',
                        'B': 'Balanced'}

        # Determine y_max
        y_max_d = 0.25
        for i_cat in i_cats:
            if y_max_d < max(rp_type_dens_dict[i_cat][e_cat].values()):
                y_max_d = max(rp_type_dens_dict[i_cat][e_cat].values())
        y_max_d = y_max_d + y_max_d / 10

        x = np.arange(len(rp_cat_labels))  # the label locations
        width = 0.35  # the width of the bars
        density_tick_labels = np.arange(0, 1, 0.25)

        fig, ax = plt.subplots(len(i_cats), figsize=(4.5, len(i_cats) * (7 / 3)))

        row_idx = 0
        for i_cat in i_cats:
            # Create bar chart
            ax[row_idx].bar(x, rp_type_dens_dict[i_cat][e_cat].values(), width, color=rp_cat_colors)
            ax[row_idx].set_title(i_cat_titles[i_cat], loc='left')
            ax[row_idx].set_xticks(x)
            ax[row_idx].set_xticklabels(rp_cat_labels)
            ax[row_idx].axhline(0.25, zorder=0, color='gray', linewidth=0.5)
            ax[row_idx].set_xlabel('Read pair type')
            ax[row_idx].set_ylabel('Density')
            ax[row_idx].set_yticks(ticks=density_tick_labels)

            # Add second axis with absolute counts and normalize to maximum density
            ax2 = ax[row_idx].twinx()
            ax2.bar(x, rp_type_freq_dict[i_cat][e_cat].values(), width, color=rp_cat_colors, edgecolor=None,
                    linewidth=3)
            ax2.set_ylabel('Read count')
            ax[row_idx].set_ylim(0, y_max_d)
            y_lim = y_max_d * sum(rp_type_freq_dict[i_cat][e_cat].values())
            ax2.set_ylim(0, y_lim)
            ax2.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

            row_idx += 1

        # Save and return plot
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    @staticmethod
    def count_configurations(d11_interaction_set: DiachromaticInteractionSet = None):

        # Initialize dictionary for counts
        conf_freq_dict = dict()
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B']:
            conf_freq_dict[i_cat] = dict()
            for e_cat in ['NN', 'EE', 'NE', 'EN', 'ALL']:
                conf_freq_dict[i_cat][e_cat] = dict()
                for i_conf in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                    conf_freq_dict[i_cat][e_cat][i_conf] = 0

        # Get absolute frequencies
        for d11_inter in d11_interaction_set.interaction_list:
            i_cat = d11_inter.get_category()
            e_cat = d11_inter.enrichment_status_tag_pair
            ht_tag = d11_inter.get_ht_tag()
            conf_freq_dict[i_cat][e_cat][ht_tag] += 1
            conf_freq_dict['ALL'][e_cat][ht_tag] += 1
            conf_freq_dict[i_cat]['ALL'][ht_tag] += 1
            conf_freq_dict['ALL']['ALL'][ht_tag] += 1
            if i_cat == 'DIX' or i_cat == 'DI':
                conf_freq_dict['U'][e_cat][ht_tag] += 1
                conf_freq_dict['U']['ALL'][ht_tag] += 1
            if i_cat == 'UIR' or i_cat == 'UI':
                conf_freq_dict['B'][e_cat][ht_tag] += 1
                conf_freq_dict['B']['ALL'][ht_tag] += 1

        # Fill second dictionary with relative frequencies
        conf_dens_dict = copy.deepcopy(conf_freq_dict)
        for i_cat in ['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B']:
            for e_cat in ['NN', 'EE', 'NE', 'EN', 'ALL']:
                i_total = sum(conf_freq_dict[i_cat][e_cat].values())
                if 0 < i_total:
                    for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                        conf_dens_dict[i_cat][e_cat][ht_tag] = conf_freq_dict[i_cat][e_cat][ht_tag] / i_total
                else:
                    for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                        conf_dens_dict[i_cat][e_cat][ht_tag] = 0.0

        return conf_freq_dict, conf_dens_dict

    @staticmethod
    def print_configuration_frequency_tables(conf_freq_dict=None,
                                             conf_dens_dict=None,
                                             i_cats=['DIX', 'DI', 'UIR', 'UI', 'ALL', 'U', 'B'],
                                             e_cats=['NN', 'EE', 'NE', 'EN', 'ALL']):
        for i_cat in i_cats:
            print(i_cat)
            for e_cat in e_cats:
                print('\t\t' + e_cat, end='')
            print()
            for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                print('\t' + ht_tag + ': ', end='')
                for e_cat in e_cats:
                    print('\t' + "{:,}".format(conf_freq_dict[i_cat][e_cat][ht_tag]), end='')
                    print(' (' + "{:.2f}".format(conf_dens_dict[i_cat][e_cat][ht_tag]) + ')', end='')
                print()

    @staticmethod
    def create_configuration_frequency_bar_charts(
            conf_freq_dict=None,
            e_cat_1='ALL',
            i_cats=['DIX', 'DI', 'UIR', 'UI', 'ALL'],
            pdf_file_name='configuration_frequency_bar_charts.pdf'):
        """
        This function creates one bar chart for each interaction category.
        For each HT configuration tag there is one bar that represent the frequencies of the configurations
        in the enrichment category passed.
        """

        # Define labels and colors
        config_colors = ['grey', 'grey', 'grey', 'grey', 'pink', 'red', 'lime', 'magenta', 'blue', 'turquoise']
        i_cat_titles = {'DIX': 'Unbalanced without reference',
                        'DI': 'Unbalanced with reference',
                        'UIR': 'Balanced reference',
                        'UI': 'Unbalanced no reference',
                        'ALL': 'All',
                        'U': 'Unbalanced',
                        'B': 'Balanced'}
        xy_label_font_size = 11.2

        # Fill second dictionary with recalcultated realtive frequencies (required for second y-axis)
        ht_tag_dens_dict = copy.deepcopy(conf_freq_dict)
        for i_cat in i_cats:
            i_total = sum(conf_freq_dict[i_cat][e_cat_1].values())
            if 0 < i_total:
                for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = conf_freq_dict[i_cat][e_cat_1][ht_tag] / i_total
            else:
                for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = 0.0

        # Determine maximal densities over all interactions categories
        y_max_d = 0.00
        for i_cat in i_cats:
            if y_max_d < max(ht_tag_dens_dict[i_cat][e_cat_1].values()):
                y_max_d = max(ht_tag_dens_dict[i_cat][e_cat_1].values())
        y_max_d = y_max_d + y_max_d / 10

        # Create a figure with plots for all interaction categories
        fig, ax = plt.subplots(len(i_cats), figsize=(5, len(i_cats) * (12 / 5)))

        x_labels = ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']
        x = np.arange(len(x_labels))  # the label locations
        width = 0.35 * 2  # the width of the bars

        row_idx = 0  # Row in the plot grid
        for i_cat in i_cats:
            # The recalculated densities for the combined enrichment categories
            e1_dens = ht_tag_dens_dict[i_cat][e_cat_1].values()

            # The absolute frequencies that were passed to the function
            e1_freq = conf_freq_dict[i_cat][e_cat_1].values()

            # Create bars for densities
            ax[row_idx].bar(x, e1_dens, width, label=e_cat_1, color=config_colors, edgecolor=config_colors,
                            linewidth=1.5)

            # Add some text for labels, title and custom x-tick labels, etc.
            ax[row_idx].set_xlabel('Configuration', fontsize=xy_label_font_size)
            ax[row_idx].set_ylabel('Density', fontsize=xy_label_font_size)
            ax[row_idx].set_title(i_cat_titles[i_cat] + ' (' + i_cat + ')', loc='left')
            ax[row_idx].set_xticks(x)
            ax[row_idx].set_xticklabels(x_labels)

            # Add second axis with absolute frequencies
            ax2 = ax[row_idx].twinx()
            ax2.bar(x - 0.45 * width, e1_freq, width, color='black', edgecolor='black', linewidth=3,
                    alpha=0)  # Set to alpha=1 to control scaling of y-axes
            ax2.bar(x + 0.45 * width, e1_freq, width, color='black', edgecolor='black', linewidth=3, alpha=0)
            ax2.set_ylabel('Interaction count', fontsize=xy_label_font_size)

            # Scale second y-axis to maximum density
            ax[row_idx].set_ylim(0, y_max_d)
            y_max_f = y_max_d * (sum(e1_freq))
            ax2.set_ylim(0, y_max_f)
            ax2.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

            # Increase tick label size
            ax[row_idx].tick_params(axis='both', which='major', labelsize=11)
            ax2.tick_params(axis='both', which='major', labelsize=11)

            # Go to the next row of the plot grid
            row_idx += 1

        # Save and return plot
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig

    @staticmethod
    def create_configuration_frequency_bar_charts_2(
            conf_freq_dict=None,
            i_cats=['DIX', 'DI', 'UIR', 'UI', 'ALL'],
            e_cat_1='NE',
            e_cat_2='EN',
            e_cat_1_color='darkred',
            e_cat_2_color='darkblue',
            pdf_file_name='ht_tag_barplot_for_two_e_cats.pdf'):
        """
        This function creates a grouped bar chart for each interaction category.
        For each HT configuration tag there are two bars that represent the frequencies of configurations
        in the two enrichment categories passed.
        The sum of the 20 relative frequencies (2 enrichment states and 10 configurations) is 1
        for each interaction category.
        """

        # Define colors and labels
        config_colors = ['grey', 'grey', 'grey', 'grey', 'pink', 'red', 'lime', 'magenta', 'blue', 'turquoise']
        i_cat_titles = {'DIX': 'Unbalanced without reference',
                        'DI': 'Unbalanced with reference',
                        'UIR': 'Balanced reference',
                        'UI': 'Unbalanced no reference',
                        'ALL': 'All',
                        'U': 'Unbalanced',
                        'B': 'Balanced'}
        xy_label_font_size = 11.2

        # Fill second dictionary with recalculated relative frequencies for e_cat_1 (NE) and e_cat_2 (EN) combined
        # (required for second y-axis)
        ht_tag_dens_dict = copy.deepcopy(conf_freq_dict)
        for i_cat in i_cats:
            i_total = sum(conf_freq_dict[i_cat][e_cat_1].values()) + sum(conf_freq_dict[i_cat][e_cat_2].values())
            if 0 < i_total:
                for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = conf_freq_dict[i_cat][e_cat_1][ht_tag] / i_total
                    ht_tag_dens_dict[i_cat][e_cat_2][ht_tag] = conf_freq_dict[i_cat][e_cat_2][ht_tag] / i_total
            else:
                for ht_tag in ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']:
                    ht_tag_dens_dict[i_cat][e_cat_1][ht_tag] = 0.0
                    ht_tag_dens_dict[i_cat][e_cat_2][ht_tag] = 0.0

        # Determine maximal densities over all interactions categories and e_cat_1 (NE) and e_cat_2 (EN)
        y_max_d = 0.00
        for i_cat in i_cats:
            for e_cat in [e_cat_1, e_cat_2]:
                if y_max_d < max(ht_tag_dens_dict[i_cat][e_cat].values()):
                    y_max_d = max(ht_tag_dens_dict[i_cat][e_cat].values())
        y_max_d = y_max_d + y_max_d / 10

        # Create a figure with plots for all interaction categories
        fig, ax = plt.subplots(len(i_cats), figsize=(5, len(i_cats) * (12 / 5)))

        x_labels = ['0X', '1X', '2X', '3X', '01', '02', '03', '12', '13', '23']
        x = np.arange(len(x_labels))  # the label locations
        width = 0.35  # the width of the bars

        row_idx = 0  # Row in the plot grid
        for i_cat in i_cats:
            # The recalculated densities for the two combined enrichment categories
            e1_dens = ht_tag_dens_dict[i_cat][e_cat_1].values()
            e2_dens = ht_tag_dens_dict[i_cat][e_cat_2].values()

            # The absolute frequencies that were passed to the function
            e1_freq = conf_freq_dict[i_cat][e_cat_1].values()
            e2_freq = conf_freq_dict[i_cat][e_cat_2].values()

            # Create bars for densities
            ax[row_idx].bar(x - 0.5 * width, e1_dens, width, label=e_cat_1, color=e_cat_1_color,
                            edgecolor=e_cat_1_color, linewidth=1.5)
            ax[row_idx].bar(x + 0.5 * width, e2_dens, width, label=e_cat_2, color=e_cat_2_color,
                            edgecolor=e_cat_2_color, linewidth=1.5)

            # Add some text for labels, title and custom x-tick labels, etc.
            ax[row_idx].set_xlabel('Configuration', fontsize=xy_label_font_size)
            ax[row_idx].set_ylabel('Density', fontsize=xy_label_font_size)
            ax[row_idx].set_title(i_cat_titles[i_cat] + ' (' + i_cat + ')', loc='left')
            ax[row_idx].set_xticks(x)
            ax[row_idx].set_xticklabels(x_labels)
            e1_patch = mpatches.Patch(color=e_cat_1_color, label='5\'')
            e2_patch = mpatches.Patch(color=e_cat_2_color, label='3\'')
            ax[row_idx].legend(handles=[e1_patch, e2_patch], fontsize=9.5, loc='upper left')

            # Add second axis with absolute frequencies
            ax2 = ax[row_idx].twinx()
            ax2.bar(x - 0.9 * width, e1_freq, width, color='red', edgecolor='red', linewidth=1.5,
                    alpha=0)  # Set to alpha=1 to control scaling of y-axes
            ax2.bar(x + 0.9 * width, e2_freq, width, color='red', edgecolor='red', linewidth=1.5, alpha=0)
            ax2.set_ylabel('Interaction count', fontsize=xy_label_font_size)

            # Scale second y-axis to maximum density
            ax[row_idx].set_ylim(0, y_max_d)
            y_max_f = y_max_d * (sum(e1_freq) + sum(
                e2_freq))  # This is possible, because we are using the combined densities of NE and EN
            ax2.set_ylim(0, y_max_f)
            ax2.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

            # Increase tick label sizes
            ax[row_idx].tick_params(axis='both', which='major', labelsize=11)
            ax2.tick_params(axis='both', which='major', labelsize=11)

            # Go to the next row of the plot grid
            row_idx += 1

        # Save and return plot
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig
