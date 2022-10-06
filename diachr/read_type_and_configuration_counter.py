import copy
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
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
            ax2.set_ylabel('Interaction count')
            ax[row_idx].set_ylim(0, y_max_d)
            y_lim = y_max_d * sum(rp_type_freq_dict[i_cat][e_cat].values())
            ax2.set_ylim(0, y_lim)
            ax2.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

            row_idx += 1

        # Save and return plot
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig
