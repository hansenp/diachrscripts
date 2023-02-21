import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns


class UnbaitedFragmentAnalysis:
    """
    This class contains functions required for the analysis of unbaited fragments performed in the Jupyter notebook
    'jupyter_notebooks/analysis/unbaited_fragment_analysis.ipynb'.
    """

    def __init__(self):

        self.foo = 1

    @staticmethod
    def get_summary_report(b_list=None, u_list=None, test='MWU'):
        print("Fragments involved only in balanced interactions: " + '{:,}'.format(len(b_list)))
        print("Fragments involved only in unbalanced interactions: " + '{:,}'.format(len(u_list)))
        b_mdn = np.median(b_list)
        u_mdn = np.median(u_list)
        b_mean = np.mean(b_list)
        u_mean = np.mean(u_list)
        if 1 < b_mdn:
            print("Median for fragments involved only in balanced interactions: " + '{:,}'.format(int(b_mdn)))
            print("Median for fragments involved only in unbalanced interactions: " + '{:,}'.format(int(u_mdn)))
        else:
            print("Median for fragments involved only in balanced interactions: " + '{:.2f}'.format(b_mdn))
            print("Median for fragments involved only in unbalanced interactions: " + '{:.2f}'.format(u_mdn))
        print("Test whether fragment lengths differ significantly: Balanced vs. unbalanced")
        if test == 'MWU':
            print('\t' + str(stats.mannwhitneyu(b_list, u_list, alternative='two-sided')))
            print('\t' + str(stats.mannwhitneyu(u_list, b_list, alternative='two-sided')))
        elif test == 'RKS':
            print('\t' + str(stats.ranksums(b_list, u_list, alternative='two-sided')))
            print('\t' + str(stats.ranksums(u_list, b_list, alternative='two-sided')))
        else:
            print('[ERROR] \'test\' must be \'MWU\' or \'RKS\'!')

    @staticmethod
    def create_boxplots(
            u_list=None,
            b_list=None,
            showfliers=True,
            title='TITEL',
            xlabel='XLAB_L',
            output_pdf='OUT.pdf'):
        # Set up figure
        fig, (ax1, ax2) = plt.subplots(1, 2)
        box_width = 0.5
        ax1.set_title(title, loc='left')

        # Determine range to be displayed
        xmin = min(u_list + b_list)
        xmax = max(u_list + b_list)
        padding = (xmax - xmin) / 30
        xmin = xmin - padding
        xmax = xmax + padding

        # Create boxplots
        bp1 = ax1.boxplot(
            [u_list, b_list],
            widths=(box_width, box_width),
            patch_artist=True,
            labels=['Unbalanced\nn=' + '{:,}'.format(len(u_list)),
                    'Balanced\nn=' + '{:,}'.format(len(b_list))],
            vert=False,
            showfliers=showfliers
        )
        colors = ['orange', 'lightblue']
        for patch, color in zip(bp1['boxes'], colors):
            patch.set_facecolor(color)
        ax1.set_xlabel(xlabel)

        # Highlight the area between the first and the third quantiles of the 5' distances of BDC0 in gray
        shaded_q1, shaded_q2, shaded_q3 = np.quantile(b_list, [0.25, 0.50, 0.75])
        ax1.axvspan(shaded_q1, shaded_q3, facecolor='gray', alpha=0.2)
        ax1.axvline(shaded_q2, color='lightgray', zorder=0)

        # Set limits of the x-axes
        # if showfliers is False:
        xmin_1, xmax_1 = ax1.get_xlim()
        xmin = min([xmin_1, xmax_1])
        xmax = max([xmin_1, xmax_1])
        ax1.set_xlim(xmin, xmax)
        padding = (xmax - xmin) / 30
        xmin = xmin - padding
        xmax = xmax + padding

        # Format figure and write to PDF file
        fig.set_figheight(3)
        fig.set_figwidth(8.5)  # 7.5
        fig.tight_layout()
        fig.savefig(output_pdf)

    @staticmethod
    def create_histplots(
            b_list=list(),
            u_list=list(),
            q_lim=1.00,
            n_bins=50,
            xlabel='XLAB_L',
            output_pdf='OUT.pdf'):

        fig_num = 5
        fig_width = 5.35
        fig_height = fig_num * (fig_width / 1.9)
        fig, ax = plt.subplots(nrows=fig_num, ncols=1, figsize=(fig_width, fig_height))

        # Compute median
        u_mdn = np.median(u_list)
        b_mdn = np.median(b_list)
        u_mad = stats.median_abs_deviation(u_list)
        b_mad = stats.median_abs_deviation(b_list)

        # Get bins
        q_max_b = np.quantile(b_list, q_lim)
        q_max_u = np.quantile(u_list, q_lim)
        x_max = max(b_list + u_list)
        if 1 < x_max:
            xq_max = max(q_max_b, q_max_u)
            bin_width = int(xq_max / n_bins)
            bins = range(0, x_max + bin_width, bin_width)
        else:
            n_bins = 25
            bin_width = 1 / n_bins
            bins = np.arange(0, 1 + bin_width, bin_width)
            xq_max = 1

        # Balanced fragment counts
        sns.histplot(ax=ax[0],
                     x=b_list,
                     stat='count',
                     bins=bins,
                     color='lightblue',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     alpha=1,
                     label='n: ' + '{:,}'.format(len(b_list)))
        ax[0].set_title('Fragments involved in balanced interactions', loc='left')
        ax[0].set_xlabel(xlabel)
        ax[0].axvline(b_mdn, linestyle='--', linewidth=1.0, color='blue')
        ax[0].axvspan(b_mdn, b_mdn + b_mad, color='green', alpha=0.25, zorder=0)
        ax[0].legend(fontsize=9)
        ax[0].set_xlim(0, xq_max)

        # Unbalanced fragment counts
        sns.histplot(ax=ax[1],
                     x=u_list,
                     stat='count',
                     bins=bins,
                     color='orange',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     alpha=1,
                     label='n: ' + '{:,}'.format(len(u_list)))
        ax[1].set_title('Fragments involved in unbalanced interactions', loc='left')
        ax[1].set_xlabel(xlabel)
        ax[1].axvline(u_mdn, linestyle='--', linewidth=1.0, color='blue')
        ax[1].axvspan(u_mdn, u_mdn + u_mad, color='green', alpha=0.25, zorder=0)
        ax[1].legend(fontsize=9)
        ax[1].set_xlim(0, xq_max)

        # Balanced fragment densities
        sns.histplot(ax=ax[2],
                     x=b_list,
                     stat='probability',
                     bins=bins,
                     color='lightblue',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     alpha=1,
                     label='n: ' + '{:,}'.format(len(b_list)))
        ax[2].set_title('Fragments involved in balanced interactions', loc='left')
        ax[2].set_xlabel(xlabel)
        ax[2].set_ylabel('Density')
        ax[2].axvline(b_mdn, linestyle='--', linewidth=1.0, color='blue')
        ax[2].axvspan(b_mdn, b_mdn + b_mad, color='green', alpha=0.25, zorder=0)
        ax[2].legend(fontsize=9)
        ax[2].set_xlim(0, xq_max)

        # Unbalanced fragment densities
        sns.histplot(ax=ax[3],
                     x=u_list,
                     stat='probability',
                     bins=bins,
                     color='orange',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     alpha=1,
                     label='n: ' + '{:,}'.format(len(u_list)))
        ax[3].set_title('Fragments involved in unbalanced interactions', loc='left')
        ax[3].set_xlabel(xlabel)
        ax[3].set_ylabel('Density')
        ax[3].axvline(u_mdn, linestyle='--', linewidth=1.0, color='blue')
        ax[3].axvspan(u_mdn, u_mdn + u_mad, color='green', alpha=0.25, zorder=0)
        ax[3].legend(fontsize=9)
        ax[3].set_xlim(0, xq_max)

        # Density comparison for balanced and unbalanced
        sns.histplot(ax=ax[4],
                     x=u_list,
                     stat='probability',
                     bins=bins,
                     color='orange',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     label='n: ' + '{:,}'.format(len(u_list)))
        sns.histplot(ax=ax[4],
                     x=b_list,
                     stat='probability',
                     bins=bins,
                     color='lightblue',
                     edgecolor="dimgray",
                     linewidth=0.5,
                     label='n: ' + '{:,}'.format(len(b_list)))
        ax[4].set_title('Density comparison for balanced and unbalanced', loc='left')
        ax[4].set_xlabel(xlabel)
        ax[4].set_ylabel('Density')
        ax[4].legend(fontsize=9)
        ax[4].set_xlim(0, xq_max)

        fig.tight_layout(pad=1.25)
        fig.savefig(output_pdf)

