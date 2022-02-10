from collections import defaultdict
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.diachromatic_interaction import DiachromaticInteraction
from typing import List


class CHCTadViz:
    """
    This class coordinates the visualization of TADs.
    """

    def __init__(self, i_file: str = None, t_file: str = None, b_file: str = None, verbose: bool = False) -> None:
        if not isinstance(i_file, str) or not os.path.isfile(i_file):
            raise ValueError("Invalid path {}".format(i_file))

        # Set plot colors
        self._i_cat_colors = {
            'DIX': 'orangered',
            'DI': 'orange',
            'UIR': 'green',
            'UI': 'gray'
        }
        self._ht_tag_colors = {
            '01': 'pink',
            '02': 'red',
            '03': 'lime',
            '12': 'magenta',
            '13': 'blue',
            '23': 'turquoise'
        }

        # Global variables for converting genomic coordinates to plot coordinates and vice versa
        self._figure_size = 10
        self._min = None
        self._max = None
        self._span = None
        self._factor = None

        if verbose:
            print('[INFO] Initializing DirectedTadVisualizer object ...')
            print('\t[INFO] Reading interactions and group them by chromosomes ...')

        self._d_inter_by_chrom_dict = dict()
        d_inter_set = DiachromaticInteractionSet(rpc_rule='ht')
        d_inter_set.parse_file(i_file=i_file, verbose=False)
        for d_inter in d_inter_set.interaction_list:
            if d_inter.chrA == d_inter.chrB:
                if d_inter.chrA not in self._d_inter_by_chrom_dict:
                    self._d_inter_by_chrom_dict[d_inter.chrA] = [d_inter]
                else:
                    self._d_inter_by_chrom_dict[d_inter.chrA].append(d_inter)

        if verbose:
            for key, d_inter_list in self._d_inter_by_chrom_dict.items():
                print('\t\t' + key + ': ' + '{:,}'.format(len(d_inter_list)) + ' interactions')
            print('\t[INFO] ... done.')

        if t_file is None:
            if verbose:
                print('[INFO] ... done.')
            return

        if verbose:
            print('\t[INFO] Reading TAD regions and group them by chromosomes ...')

        self._tads_by_chrom_dict = dict()
        with open(t_file, 'rt') as fp:
            next(fp)
            for line in fp:
                c, s, e = line.rstrip().split('\t')

                if c not in self._tads_by_chrom_dict:
                    self._tads_by_chrom_dict[c] = [(int(s), int(e))]
                else:
                    self._tads_by_chrom_dict[c].append((int(s), int(e)))

        if verbose:
            for key, tad_list in self._tads_by_chrom_dict.items():
                print('\t\t' + key + ': ' + '{:,}'.format(len(tad_list)) + ' TADs')
            print('\t[INFO] ... done.')

        if b_file is None:
            if verbose:
                print('[INFO] ... done.')
            return

        if verbose:
            print('\t[INFO] Reading baited digest regions and group them by chromosomes ...')

        self._baits_by_chrom_dict = dict()
        with open(b_file, 'rt') as fp:
            next(fp)
            for line in fp:
                c, s, e = line.rstrip().split('\t')

                if c not in self._baits_by_chrom_dict:
                    self._baits_by_chrom_dict[c] = [(int(s), int(e))]
                else:
                    self._baits_by_chrom_dict[c].append((int(s), int(e)))

        if verbose:
            for key, bait_list in self._baits_by_chrom_dict.items():
                print('\t\t' + key + ': ' + '{:,}'.format(len(bait_list)) + ' Baits')
            print('\t[INFO] ... done.')
            print('[INFO] ... done.')

    def extract_interactions(self, chrom: str, begin: int, end: int):
        """
        Given genomic coordinates, this function returns a list of Diachromatic interactions that are completely
        contained in the region to be visualized.
        """

        # Check arguments
        if not isinstance(chrom, str):
            raise ValueError("chrom must be a string")
        if not isinstance(begin, int) or not isinstance(end, int):
            raise ValueError("begin and end must be integers")

        # If there are no interactions for this chromosome, return empty list
        if chrom not in self._d_inter_by_chrom_dict:
            return []

        # Extract interactions
        inter_list = []
        for d_inter in self._d_inter_by_chrom_dict[chrom]:

            # Only keep interactions completely within the region we want to visualize
            if d_inter.chrA == chrom and begin < d_inter.fromA and d_inter.toB < end:
                inter_list.append(d_inter)

        return inter_list

    def pos_to_coordinate(self, pos):
        """
        Convert genomic position to plot coordinate.
        """
        return (pos - self._min) * self._factor

    def coordinate_to_pos(self, coordinate):
        """
        Convert plot coordinate to genomic position.
        """
        return coordinate / self._factor + self._min

    def tad_to_grey_triangle(self, start: int, end: int) -> PolygonPatch:
        begin = self.pos_to_coordinate(start)
        end = self.pos_to_coordinate(end)
        TANGENT_45 = math.tan(math.pi / 4)
        # bc is the point where lines from b and c meet
        midpoint_x = begin + 0.5 * (end - begin)
        # bc is the length of the segment from b to bc_x
        bc = midpoint_x - begin
        # tan = opposite/adjacent
        # opposite = tan(45)*adjacent
        bc_y = TANGENT_45 * bc
        poly = Polygon([(begin, 0), (midpoint_x, bc), (end, 0)])
        color = 'gray'
        return PolygonPatch(polygon=poly, color=color, alpha=0.25)

    def black_triangle_left(self, start: int, end: int):
        begin = self.pos_to_coordinate(start)
        end = self.pos_to_coordinate(end)
        midpoint_x = begin + 0.5 * (end - begin)
        bc = midpoint_x - begin
        poly = Polygon([(begin, 0), (begin, bc), (midpoint_x, bc)])
        color = 'black'
        return PolygonPatch(polygon=poly, color=color)

    def black_triangle_right(self, start: int, end: int):
        begin = self.pos_to_coordinate(start)
        end = self.pos_to_coordinate(end)
        midpoint_x = begin + 0.5 * (end - begin)
        bc = midpoint_x - begin
        poly = Polygon([(midpoint_x, bc), (end, bc), (end, 0)])
        color = 'black'
        return PolygonPatch(polygon=poly, color=color)

    def interaction_to_polygon(self, d_inter, pp_color, pp_alpha, d_radius: int = 0):
        """
        Creates a PolygonPatch for a given interaction.
        """
        if not isinstance(d_inter, DiachromaticInteraction):
            raise ValueError("Not a DiachromaticInteraction")

        # a,b,c,d are coordinates of the two digests on the X axis
        if 0 < d_radius:
            center_1 = d_inter.fromA + (d_inter.toA - d_inter.fromA) / 2
            a = self.pos_to_coordinate(center_1 - d_radius)
            b = self.pos_to_coordinate(center_1 + d_radius)
            center_2 = d_inter.fromB + (d_inter.toB - d_inter.fromB) / 2
            c = self.pos_to_coordinate(center_2 - d_radius)
            d = self.pos_to_coordinate(center_2 + d_radius)
            #print(str(center_1) + '\t' + str(a) + '\t' + str(b) + '\t' + str(c) + '\t' + str(d))
        else:
            a = self.pos_to_coordinate(d_inter.fromA)
            b = self.pos_to_coordinate(d_inter.toA)
            c = self.pos_to_coordinate(d_inter.fromB)
            d = self.pos_to_coordinate(d_inter.toB)

        TANGENT_45 = math.tan(math.pi / 4)  # The lines go up at 45 degree angles
        # bc is the point where lines from b and c meet
        bc_x = b + 0.5 * (c - b)
        # bc is the length of the segment from b to bc_x
        bc = bc_x - b
        # tan = opposite/adjacent
        # opposite = tan(45)*adjacent
        bc_y = TANGENT_45 * bc
        # ad is the point where lines from a and d meet
        ad_x = a + 0.5 * (d - a)
        ad = ad_x - a
        ad_y = TANGENT_45 * ad
        # ac is the point where lines from a and c meet
        ac_x = a + 0.5 * (c - a)
        ac = ac_x - a
        ac_y = TANGENT_45 * ac
        # bd is the point where lines from b and d meet
        bd_x = b + 0.5 * (d - b)
        bd = bd_x - b
        bd_y = TANGENT_45 * bd
        poly = Polygon([(bc_x, bc_y), (bd_x, bd_y), (ad_x, ad_y), (ac_x, ac_y)])
        return PolygonPatch(polygon=poly, color=pp_color, alpha=pp_alpha)#, linewidth=0)

    def create_visualization(self,
                             chrom: str,
                             begin: int,
                             end: int,
                             inter_cat_list: List = ['DIX', 'DI', 'UIR', 'UI'],
                             enr_cat_list: List = ['NE', 'EN', 'EE', 'NN'],
                             ht_tag_list: List = ['01', '02', '03', '12', '13', '23'],
                             color_i_cats: bool = True,
                             d_radius: int = 0,
                             plot_title: str = 'TadViz plot',
                             pdf_file_name: str = None,
                             verbose: bool = True):
        """
        This function creates the CHCTadViz plot for given genomic coordinates.
        """
        if verbose:
            print("[INFO] Creating visualization ...")

        # Check arguments
        if not isinstance(inter_cat_list, list):
            print("[ERROR] Interaction categories must be passed as a list!")
            return
        if not isinstance(enr_cat_list, list):
            print("[ERROR] Enrichment states must be passed as a list!")
            return
        if not isinstance(ht_tag_list, list):
            print("[ERROR] HT tags must be passed as a list!")
            return

        figure_height = 6
        figure_width = 2 * figure_height

        # Get list of all interactions completely within the region to be visualized
        inter_list = self.extract_interactions(chrom=chrom,
                                               begin=begin,
                                               end=end)
        if verbose:
            print("\t[INFO] Extracted a total number of " + '{:,}'.format(len(inter_list)) + " interactions in range:")
            print('\t\t' + chrom + ':' + str(begin) + '-' + str(end))

        # If there are no interaction within this region, do nothing
        if len(inter_list) == 0:
            print("[ERROR] No interactions to be visualized!")
            return

        # Get TAD regions for this chromosome
        tads = self._tads_by_chrom_dict[chrom]

        # Get bait regions for this chromosome
        baits = self._baits_by_chrom_dict[chrom]

        # Get list of interactions to be visualized and collect read pair counts
        d_inter_list = []
        rp_total_list = []
        max_i_dist = 0
        for d_inter in inter_list:
            i_cat = d_inter.get_category()
            e_cat = d_inter.enrichment_status_tag_pair
            ht_tag = d_inter.get_ht_tag()
            if i_cat in inter_cat_list and e_cat in enr_cat_list and ht_tag in ht_tag_list:
                d_inter_list.append(d_inter)
                if max_i_dist < d_inter.i_dist:
                    max_i_dist = d_inter.i_dist
                rp_total_list.append(d_inter.rp_total)  # (for all interactions or visualized interactions only?)
                # rp_total_list.append(d_inter._log10_pval)

        # If there are no interaction within this region, do nothing
        if len(d_inter_list) == 0:
            print("[ERROR] After filtering, there are no interactions left to be visualized!")
            return

        if verbose:
            print("\t[INFO] Filter for interactions to be visualized:")
            print("\t\t[INFO] Interaction categories:")
            print("\t\t\t" + str(inter_cat_list))
            print("\t\t[INFO] Enrichment status:")
            print("\t\t\t" + str(enr_cat_list))
            print("\t\t[INFO] HT tag:")
            print("\t\t\t" + str(ht_tag_list))
            print("\t\t[INFO] " + "Number interactions:")
            print("\t\t\t" + '{:,}'.format(len(d_inter_list)))

        # Determine read pair counts of interactions to be visualized for quantiles 0.1, ..., 1.0
        min_q = 0.1
        q_num = 10
        quantile_range = np.arange(0.1, 1.1, 0.1)
        quantile_values = np.quantile(rp_total_list, quantile_range)

        if verbose:
            print("\t[INFO] Read pair count quantiles of interactions to be visualized:")
            print("\t\t" + str(str(quantile_range)))
            print("\t\t" + str(quantile_values))

        # Set scaling factor for converting genomic coordinates to plot coordinates and vice versa
        self._min = begin
        self._max = end
        self._span = end - begin
        self._factor = self._figure_size / self._span

        # Create figure
        fig, ax = plt.subplots(1, figsize=(figure_width, figure_height))
        ax.set_title(plot_title, loc='left', fontsize='x-large')
        ax.set_xlabel('Genomic coordinate', labelpad=12, fontsize=14.5)
        ax.set_ylabel('Interaction distance', labelpad=12, fontsize=14.5)

        # Set limits of x and y axes
        xrange = [self.pos_to_coordinate(begin), self.pos_to_coordinate(end)]
        yrange = [0, 0.5 * (self.pos_to_coordinate(end) - self.pos_to_coordinate(begin))]
        # Uncomment the following line in order to set the upper ylim to the maximum interaction distance
        # yrange = [0, 0.5 * (self.pos_to_coordinate(2*max_i_dist) - self.pos_to_coordinate(max_i_dist))]
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)

        # Plot the TADs
        for (left_boundary, right_boundary) in tads:
            polypatch = self.tad_to_grey_triangle(start=left_boundary, end=right_boundary)
            ax.add_patch(polypatch)

        # Plot the baits
        for (sta_pos, end_pos) in baits:
            #begin = self.pos_to_coordinate(sta_pos)
            #end = self.pos_to_coordinate(end_pos)
            x = self.pos_to_coordinate(sta_pos) + (self.pos_to_coordinate(end_pos) - self.pos_to_coordinate(sta_pos)) / 2
            ax.axvline(x, color='gray', linewidth=0.5, zorder=0, alpha=0.5)

        # Plot the interaction polygons
        for d_inter in d_inter_list:

            # Choose color depending on the interaction category or ht_tag
            if color_i_cats:
                pp_color = self._i_cat_colors[d_inter.get_category()]
            else:
                pp_color = self._ht_tag_colors[d_inter.get_ht_tag()]

            # Determine transparency depending on read pair count
            pp_alpha = min_q
            for i in range(0, q_num):
                if quantile_values[i] <= d_inter.rp_total:
                    pp_alpha = quantile_range[i]
            # if pp_alpha >= 0.8:
            # print(str(d_inter.rp_total) + '\t' + str(pp_alpha))

            # Get and plot interaction polypatch
            polypatch = self.interaction_to_polygon(d_inter=d_inter,
                                                    pp_color=pp_color,
                                                    pp_alpha=pp_alpha,
                                                    d_radius=d_radius)
            ax.add_patch(polypatch)

        # Plot black triangles
        polypatch = self.black_triangle_left(start=begin, end=end)
        ax.add_patch(polypatch)
        polypatch = self.black_triangle_right(start=begin, end=end)
        ax.add_patch(polypatch)

        # Add genomic coordinate labels to x-axis
        ax.set_xticks(ax.get_xticks())
        xtick_labels = ['{:,}'.format(round(self.coordinate_to_pos(x))) for x in ax.get_xticks()]
        ax.set_xticklabels(xtick_labels, fontsize=11.5)

        # Add interaction distance labels to y-axis (correct?)
        ax.set_yticks(ax.get_yticks())
        ytick_labels = ['{:,}'.format(2 * round(self.coordinate_to_pos(y) - begin)) for y in ax.get_yticks()]
        ytick_labels[0] = ''
        ax.set_yticklabels(ytick_labels, fontsize=11.5)
        ax.set_aspect(1)

        if verbose:
            print("[INFO] ... done.")

        # Add annotations
        label_text = chrom + ':' + str(begin) + '-' + str(end) + '\n'
        label_text += '# Interactions: ' + '{:,}'.format(len(d_inter_list)) + '\n'
        label_text += 'I cats: ' + str(inter_cat_list) + '\n'
        label_text += 'E states: ' + str(enr_cat_list) + '\n'
        label_text += 'HT tags: ' + str(ht_tag_list)
        label_x_pos = xrange[0] + (xrange[1] - xrange[0])/60
        label_y_pos = yrange[1] - yrange[1]/4.7
        ax.annotate(label_text, xy=(label_x_pos, label_y_pos), color='white', size=12)

        # Save and return figure
        fig.tight_layout()
        if pdf_file_name is not None:
            fig.savefig(pdf_file_name)
        return fig
