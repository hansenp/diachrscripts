import os
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from re import sub


class BaitAnalysis:
    """
    This class contains functions required for the bait analysis performed in the Jupyter notebook:

        jupyter_notebooks/analysis/bait_analysis.ipynb
    """

    def __init__(self):

        # Dictionary with baited fragments and associated baits
        self.FRAGS_WITH_BAITS_dict = dict()

        # Dictionary with bait sequences
        self.BAIT_SEQ_dict = dict()

    def init_fragments_with_baits_dict(
            self,
            bait_bed_file: str = None,
            baited_fragment_bed_file: str = None,
            bfc0_bed_file: str = None,
            bfc1_bed_file: str = None,
            bfc2_bed_file: str = None,
            overlap_fraction: float = 1.00,
            verbose: bool = False):
        """
        Creates a dictionary from which, given the coordinates of a baited fragment, the coordinates of the associated
        baits and the BF class can be retrieved.
        """

        print('[INFO] Init \'FRAGS_WITH_BAITS_dict\' ...')

        self.FRAGS_WITH_BAITS_dict = dict()

        # Assign baits to fragments
        # -------------------------

        print('\t[INFO] Assign baits to their fragments ...')

        # Use PyBedTools to assign baits to their fragments
        # The option 'F=1.00' causes an overlap to be reported only if a bait is completely contained in the digest.
        bait_bt_file = pybedtools.BedTool(bait_bed_file)
        baited_fragment_bt_file = pybedtools.BedTool(baited_fragment_bed_file)
        fragments_with_bait_coords = baited_fragment_bt_file.intersect(bait_bt_file, wa=True, wb=True, F=overlap_fraction)

        # Check bait to fragment assignment
        print('\t\t[INFO] Total number of baits: ' + '{:,}'.format(len(bait_bt_file)))
        print('\t\t[INFO] Assigned baits: ' + '{:,}'.format(len(fragments_with_bait_coords)))
        print('\t\t[INFO] Unassigned baits: ' + '{:,}'.format(len(bait_bt_file) - len(fragments_with_bait_coords)))

        # Find fragments to which no baits was assigned
        frag_set_1 = set()
        for d in baited_fragment_bt_file:
            frag_set_1.add(d.chrom + '\t' + str(d.start) + '\t' + str(d.end))
        frag_set_2 = set()
        for d in fragments_with_bait_coords:
            frag_set_2.add(d.chrom + '\t' + str(d.start) + '\t' + str(d.end))
        n_frags_without_bait = len(frag_set_1) - len(frag_set_2)
        print('\t\t[INFO] Total number of baited fragments: ' + '{:,}'.format(len(frag_set_1)))
        print('\t\t[INFO] Number of fragments with at least one bait: ' + '{:,}'.format(len(frag_set_2)))
        print('\t\t[INFO] Number of fragments without bait: ' + '{:,}'.format(n_frags_without_bait))
        if verbose and 0 < n_frags_without_bait:
            print('\t\t\t[INFO] Fragments without bait:')
            frag_set_3 = frag_set_1 - frag_set_2
            for d in frag_set_3:
                print('\t\t\t\t' + d)

        # Assign fragments to their BF classes and baits to 5' or 3' ends
        # ---------------------------------------------------------------

        print('\t[INFO] Assign fragments to their BF classes and baits to 5\' or 3\' ends ...')

        # Create three sets of coordinates of BFC fragments for quick lookups
        bfc0_bt_file = pybedtools.BedTool(bfc0_bed_file)
        bfc1_bt_file = pybedtools.BedTool(bfc1_bed_file)
        bfc2_bt_file = pybedtools.BedTool(bfc2_bed_file)
        bfc0_set = set()
        for d in bfc0_bt_file:
            bfc0_set.add(d.chrom + '\t' + str(d.start) + '\t' + str(d.end))
        bfc1_set = set()
        for d in bfc1_bt_file:
            bfc1_set.add(d.chrom + '\t' + str(d.start) + '\t' + str(d.end))
        bfc2_set = set()
        for d in bfc2_bt_file:
            bfc2_set.add(d.chrom + '\t' + str(d.start) + '\t' + str(d.end))

        # Iterate through the file created with BedTools intersect
        for d in fragments_with_bait_coords:

            # Split line with digest and bait coordinates
            arr = str(d).split('\t')
            d_chr = arr[0]
            d_sta = int(arr[1])
            d_end = int(arr[2])
            b_chr = arr[3]
            b_sta = int(arr[4])
            b_end = int(arr[5])

            # Calculate the center positions of the digest and the overlapping bait
            d_len = d_end - d_sta + 1
            d_center_pos = d_sta + int(d_len / 2) + 1
            b_len = b_end - b_sta + 1
            b_center_pos = b_sta + int(b_len / 2) + 1

            # Get keys for fragment and overlapping bait
            d_key = d_chr + '\t' + str(d_sta) + '\t' + str(d_end)
            b_key = b_chr + '\t' + str(b_sta) + '\t' + str(b_end)

            # Create dictionary for this fragment if none exists yet
            if d_key not in self.FRAGS_WITH_BAITS_dict:
                self.FRAGS_WITH_BAITS_dict[d_key] = dict()
                self.FRAGS_WITH_BAITS_dict[d_key]['B5_COORDS'] = []
                self.FRAGS_WITH_BAITS_dict[d_key]['B3_COORDS'] = []

            # Assign the overlapping bait to the 5' or 3' end of the fragment
            if b_center_pos < d_center_pos:
                self.FRAGS_WITH_BAITS_dict[d_key]['B5_COORDS'].append(b_key)
            else:
                self.FRAGS_WITH_BAITS_dict[d_key]['B3_COORDS'].append(b_key)

            # Assign the fragment to one of the BF classes
            if d_key in bfc1_set:
                self.FRAGS_WITH_BAITS_dict[d_key]['BFC'] = 1
            elif d_key in bfc2_set:
                self.FRAGS_WITH_BAITS_dict[d_key]['BFC'] = 2
            elif d_key in bfc0_set:
                self.FRAGS_WITH_BAITS_dict[d_key]['BFC'] = 0
            else:
                self.FRAGS_WITH_BAITS_dict[d_key]['BFC'] = 'X'  # Fragments that are neither BFC0, BFC1, nor BFC2

        # Check fragment to BF class assignment
        cnt_dict = {'0': 0, '1': 0, '2': 0, 'X': 0}
        for d_key in self.FRAGS_WITH_BAITS_dict.keys():
            cnt_dict[str(self.FRAGS_WITH_BAITS_dict[d_key]['BFC'])] += 1
        for key in cnt_dict.keys():
            frag_num = cnt_dict[key]
            percentage = 100 * cnt_dict[key] / sum(cnt_dict.values())
            print('\t\t[INFO] BFC' + key + ': ' + '{:,}'.format(frag_num) + ' (' + '{:.2f}'.format(percentage) + '%)')
        print('\t\t[INFO] Total number of baited fragments: ' + '{:,}'.format(sum(cnt_dict.values())))
        n_bf012 = sum(cnt_dict.values()) - cnt_dict['X']
        print('\t\t[INFO] Number of BFC0, BFC1, and BFC2 fragments: ' + '{:,}'.format(n_bf012))

        # Remove BFCX fragments from dictionary
        print('\t[INFO] Remove BFCX fragments from dictionary ...')
        tmp_dict = dict()
        for d_key in self.FRAGS_WITH_BAITS_dict.keys():
            if self.FRAGS_WITH_BAITS_dict[d_key]['BFC'] != 'X':
                tmp_dict[d_key] = self.FRAGS_WITH_BAITS_dict[d_key]
        self.FRAGS_WITH_BAITS_dict = tmp_dict

        # Determine numbers of baited fragments with 1, 2 or more baits
        print('\t[INFO] Determine numbers of baited fragments with 1, 2 or more baits ...')

        cnt_dict = {'1': 0, '2': 0, '2<': 0}
        for d_key in self.FRAGS_WITH_BAITS_dict.keys():
            bait_num = len(
                self.FRAGS_WITH_BAITS_dict[d_key]['B5_COORDS'] + self.FRAGS_WITH_BAITS_dict[d_key]['B3_COORDS'])
            if 2 < bait_num:
                cnt_dict['2<'] += 1
            else:
                cnt_dict[str(bait_num)] += 1

        # Report counts
        for key in cnt_dict.keys():
            frag_num = cnt_dict[key]
            percentage = 100 * cnt_dict[key] / sum(cnt_dict.values())
            print(
                '\t\t[INFO] ' + key + ' baits: ' + '{:,}'.format(frag_num) + ' (' + '{:.2f}'.format(percentage) + '%)')
        print('\t\t[INFO] Total number of baited fragments: ' + '{:,}'.format(sum(cnt_dict.values())))

        print('[INFO] ... done.')
        return self.FRAGS_WITH_BAITS_dict

    @staticmethod
    def get_distances_between_baits_and_their_restriction_sites(bf_coords: str, b5_coords: str, b3_coords: str):

        # Get relevant start and end coordinates of fragment and baits
        bf_sta = int(bf_coords.split('\t')[1])
        bf_end = int(bf_coords.split('\t')[2])
        b5_sta = int(b5_coords.split('\t')[1])
        b3_end = int(b3_coords.split('\t')[2])

        # Calculate and return distances
        return b5_sta - bf_sta, bf_end - b3_end

    def split_frags_with_baits_dict_by_bait_characteristics(
            self,
            frags_with_baits_dict: dict() = None,
            max_dist: int = 0):

        # Count unilaterally, bilaterally shifted, and bilaterally unshifted baited fragments separately by BF class
        split_cnt_dict = {
            'BFC0': {'UNILATERAL': 0, 'BILATERAL_SHIFTED': 0, 'BILATERAL_UNSHIFTED': 0},
            'BFC1': {'UNILATERAL': 0, 'BILATERAL_SHIFTED': 0, 'BILATERAL_UNSHIFTED': 0},
            'BFC2': {'UNILATERAL': 0, 'BILATERAL_SHIFTED': 0, 'BILATERAL_UNSHIFTED': 0},
        }

        # Dictionaries to be returned
        unilateral_dict = dict()
        bilateral_shifted = dict()
        bilateral_unshifted = dict()

        # Go through all baited digests
        for d_key in frags_with_baits_dict.keys():

            # Get baited fragment class
            bd_class = 'BFC' + str(frags_with_baits_dict[d_key]['BFC'])

            # Unilaterally baited fragments
            if len(frags_with_baits_dict[d_key]['B5_COORDS']) != 1 or len(
                    frags_with_baits_dict[d_key]['B3_COORDS']) != 1:
                unilateral_dict[d_key] = frags_with_baits_dict[d_key]
                split_cnt_dict[bd_class]['UNILATERAL'] += 1
                continue

            # Get distances
            b5_dist, b3_dist = self.get_distances_between_baits_and_their_restriction_sites(
                bf_coords=d_key,
                b5_coords=frags_with_baits_dict[d_key]['B5_COORDS'][0],
                b3_coords=frags_with_baits_dict[d_key]['B3_COORDS'][0])

            # Distinguish between bilaterally baited fragments with shifted and unshifted baits
            if b5_dist <= max_dist and b3_dist <= max_dist:
                bilateral_unshifted[d_key] = frags_with_baits_dict[d_key]
                split_cnt_dict[bd_class]['BILATERAL_UNSHIFTED'] += 1
            else:
                bilateral_shifted[d_key] = frags_with_baits_dict[d_key]
                split_cnt_dict[bd_class]['BILATERAL_SHIFTED'] += 1

        return unilateral_dict, bilateral_shifted, bilateral_unshifted, split_cnt_dict

    @staticmethod
    def get_latex_table_a(split_cnt_dict: dict() = None):

        df_split_cnt = pd.DataFrame(data=split_cnt_dict)
        df_pfc = df_split_cnt.div(df_split_cnt.sum(axis=0), axis=1).round(4) * 100
        print('\nPercentages for columns:')
        print(df_pfc)
        print('\nPercentages for rows:')
        print(df_pfc.div(df_pfc.sum(axis=1), axis=0).round(4) * 100)
        print()
        df_cwt = df_split_cnt.copy()
        df_cwt.loc['Total'] = df_cwt.sum(numeric_only=True, axis=0)
        df_cwt.loc[:, 'Total'] = df_cwt.sum(numeric_only=True, axis=1)
        print('Counts with totals:')
        print(df_cwt)
        print('\nPercentages for columns:')
        total_col = df_cwt['Total']
        total_col.drop("Total", axis=0, inplace=True)
        df_pfc['Total'] = total_col.div(total_col.sum()).round(4) * 100
        print(df_pfc)
        print()
        df_cwt.rename(index={"UNILATERAL": "I. Unilateral", "BILATERAL_SHIFTED": "II. Shifted",
                             "BILATERAL_UNSHIFTED": "III. Unshifted"}, inplace=True)
        df_pfc.rename(index={"UNILATERAL": "I. Unilateral", "BILATERAL_SHIFTED": "II. Shifted",
                             "BILATERAL_UNSHIFTED": "III. Unshifted"}, inplace=True)
        print("\\begin{tabular}{l|rr|rr|rr|rr}")
        print(
            "\t\\multicolumn{9}{l}{{\\bf A:} Unilaterally baited fragments and bilaterally baited fragments with and "
            "without shifted baits}" + "\\" + "\\" + "\\hline")
        print("\t& {\\bf BFC0}& & {\\bf BFC1}& & {\\bf BFC2}& & {\\bf Total}" + "\\" + "\\" + "\\hline")
        for i in ['I. Unilateral', 'II. Shifted', 'III. Unshifted', 'Total']:
            row = "\t{\\bf " + i.replace("_", " ") + "}"
            for j in ['BFC0', 'BFC1', 'BFC2', 'Total']:
                row += "& " + "{:,}".format(df_cwt[j][i])
                if i != 'Total':
                    row += "& " + str(int(round(df_pfc[j][i]))) + "\\%"
                else:
                    if j != 'Total':
                        row += "& "
            row += "\\" + "\\"
            if i == 'III. Unshifted' or i == 'Total':
                row += "\\hline"
            print(row)
        print("\\end{tabular}")

    @staticmethod
    def get_latex_table_b(split_cnt_dict: dict() = None):

        df_split_cnt = pd.DataFrame(data=split_cnt_dict)
        df_cont = df_split_cnt.drop(columns=['BFC1', 'BFC2'])
        df_cont['BFC12'] = df_split_cnt['BFC1'] + df_split_cnt['BFC2']
        df_cont.loc['BILATERAL'] = df_cont.loc['BILATERAL_SHIFTED'] + df_cont.loc['BILATERAL_UNSHIFTED']
        df_cont.drop(labels=['BILATERAL_SHIFTED', 'BILATERAL_UNSHIFTED'], axis=0, inplace=True)
        df_cont.loc['Total'] = df_cont.sum(numeric_only=True, axis=0)
        df_cont['Total'] = df_cont.sum(numeric_only=True, axis=1)
        df_cont.rename(index={"UNILATERAL": "I. Unilateral", "BILATERAL": "II,III. Bilateral"}, inplace=True)
        print(df_cont)
        print()
        df_cont_p = df_cont.drop(labels=['Total'], axis=0, inplace=True)
        df_cont_p = df_cont.div(df_cont.sum(axis=0), axis=1).round(4) * 100
        df_cont.loc['Total'] = df_cont.sum(numeric_only=True, axis=0)
        print(df_cont_p)
        print()
        print("\\begin{tabular}{l|rr|rr|rr}")
        print(
            "\t\\multicolumn{6}{l}{{\\bf B:} Unilaterally vs. bilaterally baited fragments}" + "\\" + "\\" + "\\hline")
        print("\t& {\\bf BFC0}& & {\\bf BFC12}& & {\\bf Total}" + "\\" + "\\" + "\\hline")
        for i in ['I. Unilateral', 'II,III. Bilateral', 'Total']:
            row = "\t{\\bf " + i.replace("_", " ") + "}"
            for j in ['BFC0', 'BFC12', 'Total']:
                row += "& " + "{:,}".format(df_cont[j][i])
                if i != 'Total':
                    row += "& " + str(int(round(df_cont_p[j][i]))) + "\\%"
                else:
                    if j != 'Total':
                        row += "& "
            row += "\\" + "\\"
            if i == 'II,III. Bilateral' or i == 'Total':
                row += "\\hline"
            print(row)
        print("\\end{tabular}")

    @staticmethod
    def get_latex_table_c(baited_5_3_cnt_dict: dict() = None):

        df = pd.DataFrame(data=baited_5_3_cnt_dict)
        df.index = ['5\' baited', '3\' baited']
        df.loc['Total'] = df.sum(numeric_only=True, axis=0)
        df.loc[:, 'Total'] = df.sum(numeric_only=True, axis=1)
        df_t = df.transpose(copy=True)
        print("\\begin{tabular}{l|r|r|r}")
        print("\t\\multicolumn{4}{l}{{\\bf C:} Unilateral separated by 5' and 3'}" + "\\" + "\\" + "\\hline")
        print("\t& {\\bf 5' baited}& {\\bf 3' baited}& {\\bf Total}" + "\\" + "\\" + "\\hline")
        for i in ['BFC0', 'BFC1', 'BFC2', 'Total']:
            row = "\t{\\bf " + i.replace("_", " ") + "}"
            for j in ['5\' baited', '3\' baited', 'Total']:
                row += "& " + "{:,}".format(df_t[j][i])
            row += "\\" + "\\"
            if i == "BFC2" or i == "Total":
                row += "\\hline"
            print(row)
        print("\\end{tabular}")

    @staticmethod
    def get_latex_table_d(split_cnt_dict: dict() = None):

        df = pd.DataFrame(data=split_cnt_dict)
        df['BFC12'] = df['BFC1'] + df['BFC2']
        df.drop("BFC1", axis=1, inplace=True)
        df.drop("BFC2", axis=1, inplace=True)
        df.drop("UNILATERAL", axis=0, inplace=True)
        df.loc[:, 'Total'] = df.sum(numeric_only=True, axis=1)
        df.loc['Total'] = df.sum(numeric_only=True, axis=0)
        df.rename(index={"BILATERAL_SHIFTED": "II. Shifted", "BILATERAL_UNSHIFTED": "III. Unshifted"}, inplace=True)
        print(df)
        print()
        df_pct = df.copy()
        df_pct.drop("Total", axis=0, inplace=True)
        # df_pct.drop("Total", axis=1, inplace=True)
        df_pct = df_pct.div(df_pct.sum(axis=0), axis=1).round(4) * 100
        print(df_pct)
        print()
        print("\\begin{tabular}{l|rr|rr|rr}")
        print(
            "\t\\multicolumn{6}{l}{{\\bf D:} Bilaterally baited fragments: shifted vs. unshifted}" + "\\" + "\\" + "\\hline")
        print("\t& {\\bf BFC0}& & {\\bf BFC12}& & {\\bf Total}" + "\\" + "\\" + "\\hline")
        for i in ['II. Shifted', 'III. Unshifted', 'Total']:
            row = "\t{\\bf " + i.replace("_", " ") + "}"
            for j in ['BFC0', 'BFC12', 'Total']:
                row += "& " + "{:,}".format(df[j][i])
                if i != 'Total':
                    row += "& " + str(int(round(df_pct[j][i]))) + "\\%"
                else:
                    if j != 'Total':
                        row += "& "
            row += "\\" + "\\"
            if i == 'III. Unshifted' or i == 'Total':
                row += "\\hline"
            print(row)
        print("\\end{tabular}")

    def get_lists_of_bait_restriction_site_distances(self, bilateral_shifted_dict: dict() = None, bf_classes=[0, 1, 2]):

        # Lists of distances between baits and their restriction sites
        b5_dist_list = []
        b3_dist_list = []

        # Iterate through baited fragments
        for f_key in bilateral_shifted_dict.keys():

            # Catch errors that occur when calling function with unilaterally baited fragments
            if len(bilateral_shifted_dict[f_key]['B5_COORDS']) == 0 or len(
                    bilateral_shifted_dict[f_key]['B3_COORDS']) == 0:
                print('[ERROR] Fragment \'' + f_key + '\' is not bilaterally baited!')
                break

            # Filter for fragments with specified BF classes
            if bilateral_shifted_dict[f_key]['BFC'] not in bf_classes:
                continue

            # Get distances
            b5_dist, b3_dist = self.get_distances_between_baits_and_their_restriction_sites(
                bf_coords=f_key,
                b5_coords=bilateral_shifted_dict[f_key]['B5_COORDS'][0],
                b3_coords=bilateral_shifted_dict[f_key]['B3_COORDS'][0])
            b5_dist_list.append(b5_dist)
            b3_dist_list.append(b3_dist)

        return b5_dist_list, b3_dist_list

    @staticmethod
    def create_bfc_boxplot(
            bfc0_b5: list() = None,  # BFC0 5' bait
            bfc0_b3: list() = None,  # BFC0 3' bait
            bfc1_b5: list() = None,  # BFC1 5' bait
            bfc1_b3: list() = None,  # BFC1 3' bait
            bfc2_b5: list() = None,  # BFC2 5' bait
            bfc2_b3: list() = None,  # BFC2 3' bait
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
        xmin = min(bfc0_b5 + bfc0_b3 + bfc1_b5 + bfc1_b3 + bfc2_b5 + bfc2_b3)
        xmax = max(bfc0_b5 + bfc0_b3 + bfc1_b5 + bfc1_b3 + bfc2_b5 + bfc2_b3)

        # Create boxplots for 5' baits (left)
        bp1 = ax1.boxplot(
            [bfc2_b5, bfc1_b5, bfc0_b5],
            widths=(box_width, box_width, box_width),
            patch_artist=True,
            labels=['BFC2\nn=' + '{:,}'.format(len(bfc2_b5)),
                    'BFC1\nn=' + '{:,}'.format(len(bfc1_b5)),
                    'BFC0\nn=' + '{:,}'.format(len(bfc0_b5))],
            vert=False,
            showfliers=showfliers
        )
        colors = ['green', 'blue', 'gray']
        for patch, color in zip(bp1['boxes'], colors):
            patch.set_facecolor(color)
        ax1.set_xlabel(xlabel_l)

        # Highlight the area between the first and the third quantiles of the 5' distances of BFC0 in gray
        shaded_q1, shaded_q2, shaded_q3 = np.quantile(bfc0_b5, [0.25, 0.50, 0.75])
        ax1.axvspan(shaded_q1, shaded_q3, facecolor='gray', alpha=0.2)
        ax1.axvline(shaded_q2, color='lightgray', zorder=0)

        # Create boxplots for 3' baits (right)
        x = ax2.boxplot(
            [bfc2_b3, bfc1_b3, bfc0_b3],
            widths=(box_width, box_width, box_width),
            patch_artist=True,
            vert=False,
            showfliers=showfliers
        )
        colors = ['green', 'blue', 'gray']
        for patch, color in zip(x['boxes'], colors):
            patch.set_facecolor(color)
        ax2.set_xlabel(xlabel_r)

        # Hide tick labels
        [t.set_color('white') for t in ax2.yaxis.get_ticklabels()]

        # Highlight the area between the first and the third quantiles of the 5' distances of BFC0 in gray
        shaded_q1, shaded_q2, shaded_q3 = np.quantile(bfc0_b3, [0.25, 0.50, 0.75])
        ax2.axvspan(shaded_q1, shaded_q3, facecolor='gray', alpha=0.2)
        ax2.axvline(shaded_q2, color='lightgray', zorder=0)

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

    def init_bait_seq_dict(self,
                           genome_fasta_file: str = None,
                           bait_bed_file: str = None,
                           working_directory: str = '.'):
        """
        Builds a dictionary with bait sequences using the bait coordinates as keys.

        :param genome_fasta_file: FASTA file with genome sequence. Inddex (.fai) with the same name must be in the
        same place.
        :param bait_bed_file: BED file with bait coordinates
        :param working_directory: Directory for temporary file
        """
        bait_seqs = pybedtools.BedTool(bait_bed_file)
        bait_seqs.sequence(fi=genome_fasta_file, fo=working_directory + '/bait_sequences_tmp.fa', tab=True)
        with open(working_directory + '/bait_sequences_tmp.fa', 'rt') as fp:
            for line in fp:
                coords, seq = line.rstrip().split('\t')
                chrom = coords.split(':')[0]
                sta = coords.split(':')[1].split('-')[0]
                end = coords.split(':')[1].split('-')[1]
                self.BAIT_SEQ_dict[chrom + '\t' + sta + '\t' + end] = seq
        os.remove(working_directory + '/bait_sequences_tmp.fa')

    def get_bait_seq(self, bait_coords: str = None):
        return self.BAIT_SEQ_dict[bait_coords]

    @staticmethod
    def get_gc_content_seq(sequence: str = None):
        """
        Returns the GC content of a given sequence.
        """
        gc_count = sequence.count('g')
        gc_count += sequence.count('G')
        gc_count += sequence.count('c')
        gc_count += sequence.count('C')
        seq_len = len(sequence)
        gc_content = gc_count / seq_len
        return gc_content

    def get_gc_content(self, coords: str = None):
        """
        Returns the GC content of a sequence at the given coordinates.
        """
        bait_seq = self.get_bait_seq(coords)
        return self.get_gc_content_seq(bait_seq)

    @staticmethod
    def get_repeat_content_seq(sequence: str = None):
        """
        Returns the repeat content of a given sequence.
        """
        sequence = sub("[a-z]", 'R', sequence)
        r_count = sequence.count('R')
        seq_len = len(sequence)
        repeat_content = r_count / seq_len
        return repeat_content

    def get_repeat_content(self, coords: str = None):
        """
        Returns the repeat content of a sequence at the given coordinates.
        """
        bait_seq = self.get_bait_seq(coords)
        return self.get_repeat_content_seq(bait_seq)

    def get_lists_of_bait_gc_content(self, bilateral_unshifted_dict: dict() = None, bf_classes=[0, 1, 2]):

        b5_gc_list = []
        b3_gc_list = []
        for f_key in bilateral_unshifted_dict.keys():

            # Catch errors that occur when calling function with unilaterally baited fragments
            if len(bilateral_unshifted_dict[f_key]['B5_COORDS']) == 0 or len(
                    bilateral_unshifted_dict[f_key]['B3_COORDS']) == 0:
                print('[ERROR] Fragment \'' + f_key + '\' is not bilaterally baited!')
                break

            # Filter for fragments with specified BF classes
            if bilateral_unshifted_dict[f_key]['BFC'] not in bf_classes:
                continue

            # Get GC content of 5' bait
            bait_5_coord = bilateral_unshifted_dict[f_key]['B5_COORDS'][0]
            gc_content = self.get_gc_content(bait_5_coord)
            b5_gc_list.append(gc_content)

            # Get GC content of 3' bait
            bait_3_coord = bilateral_unshifted_dict[f_key]['B3_COORDS'][0]
            gc_content = self.get_gc_content(bait_3_coord)
            b3_gc_list.append(gc_content)

        return b5_gc_list, b3_gc_list

    def get_lists_of_bait_repeat_content(self, bilateral_unshifted_dict: dict() = None, bf_classes=[0, 1, 2]):

        b5_repeat_list = []
        b3_repeat_list = []
        for f_key in bilateral_unshifted_dict.keys():

            # Catch errors that occur when calling function with unilaterally baited fragments
            if len(bilateral_unshifted_dict[f_key]['B5_COORDS']) == 0 or len(
                    bilateral_unshifted_dict[f_key]['B3_COORDS']) == 0:
                print('[ERROR] Fragment \'' + f_key + '\' is not bilaterally baited!')
                break

            # Filter for fragments with specified BF classes
            if bilateral_unshifted_dict[f_key]['BFC'] not in bf_classes:
                continue

            # Get repeat content of 5' bait
            bait_5_coord = bilateral_unshifted_dict[f_key]['B5_COORDS'][0]
            repeat_content = self.get_repeat_content(bait_5_coord)
            b5_repeat_list.append(repeat_content)

            # Get repeat content of 3' bait
            bait_3_coord = bilateral_unshifted_dict[f_key]['B3_COORDS'][0]
            repeat_content = self.get_repeat_content(bait_3_coord)
            b3_repeat_list.append(repeat_content)

        return b5_repeat_list, b3_repeat_list
