import re
import matplotlib.pyplot as plt
import numpy as np


class OriginalBaitAnalysis:
    """
    This class contains functions required for the analysis of the original baits used in Mifsud et al. 2015
    and Schoenfelder et al. 2015 performed in the Jupyter notebook:

        jupyter_notebooks/analysis/original_bait_analysis.ipynb
    """

    def __init__(self):
        pass

    @staticmethod
    def get_period_lengths(seq: str = None, ptype: str = 'repeat', verbose: bool = False):
        """
        Takes a sequencence of ACGTNacgtn and returns a list with the lengths of each repeat period
        """
        seq = 'A' + seq + 'A'
        if ptype == 'repeat':
            sta_positions = [m.start() for m in re.finditer('[NACGT][nacgt]', seq)]
            end_positions = [m.start() for m in re.finditer('[acgtn][NACGT]', seq)]
        elif ptype == 'N':
            sta_positions = [m.start() for m in re.finditer('[acgtACGT]N', seq)]
            end_positions = [m.start() for m in re.finditer('N[acgtACGT]', seq)]
        else:
            print("[ERROR] Invalid 'ptype' argument! Must be 'repeat' or 'N'.")
            return
        regions = list(zip(sta_positions, end_positions))
        len_list = [end - sta for sta, end in regions]
        return len_list
    
    @staticmethod
    def get_n_content_seq(sequence: str = None):
        """
        Returns the 'N' content of a given sequence.
        """
        n_count = sequence.count('N')
        seq_len = len(sequence)
        n_content = n_count / seq_len
        return n_content


    @staticmethod
    def create_contents_boxplot(
            gc_l: list = None,  # BFC0 5' bait
            gc_r: list = None,  # BFC0 3' bait
            rep_l: list = None,  # BFC1 5' bait
            rep_r: list = None,  # BFC1 3' bait
            n_l: list = None,  # BFC2 5' bait
            n_r: list = None,  # BFC2 3' bait
            showfliers: bool = True,
            title: str = 'TITLE',
            xlabel_l: str = 'XLAB_L',
            xlabel_r: str = 'XLAB_R',
            output_pdf: str = 'OUT.pdf'):

        # Set up figure
        fig, (ax1, ax2) = plt.subplots(1, 2)
        box_width = 0.5
        ax1.set_title(title, loc='left')

        # Determine range to be displayed
        xmin = min(gc_l + gc_r + rep_l + rep_r + n_l + n_r)
        xmax = max(gc_l + gc_r + rep_l + rep_r + n_l + n_r)

        # Create boxplots for 5' baits (left)
        bp1 = ax1.boxplot(
            [n_l, rep_l, gc_l],
            widths=(box_width, box_width, box_width),
            patch_artist=True,
            tick_labels=['N', 'Repeat', 'GC'],
            vert=False,
            showfliers=showfliers
        )
        colors = ['royalblue', 'crimson', 'gold']
        for patch, color in zip(bp1['boxes'], colors):
            patch.set_facecolor(color)
        ax1.set_xlabel(xlabel_l)
        ax1.grid(True)

        # Create boxplots for 3' baits (right)
        x = ax2.boxplot(
            [n_r, rep_r, gc_r],
            widths=(box_width, box_width, box_width),
            patch_artist=True,
            vert=False,
            showfliers=showfliers
        )
        for patch, color in zip(x['boxes'], colors):
            patch.set_facecolor(color)
        ax2.set_xlabel(xlabel_r)
        ax2.grid(True)

        # Hide tick labels
        [t.set_color('white') for t in ax2.yaxis.get_ticklabels()]

        # Set limits of the x-axes
        # if showfliers is False:
        xmin_1, xmax_1 = ax1.get_xlim()
        xmin_2, xmax_2 = ax2.get_xlim()
        xmin = min([xmin_1, xmax_1, xmin_2, xmax_2])
        xmax = max([xmin_1, xmax_1, xmin_2, xmax_2])
        ax1.set_xlim(xmin, xmax)
        ax2.set_xlim(xmin, xmax)

        # Format figure and write to PDF file
        fig.set_figheight(3)
        fig.set_figwidth(8.5)  # 7.5
        fig.tight_layout()
        fig.savefig(output_pdf)