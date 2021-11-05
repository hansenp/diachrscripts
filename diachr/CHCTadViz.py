from collections import defaultdict
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.diachromatic_interaction import DiachromaticInteraction


class CHCTadViz:
    """
    This class coorrdinates the visualization of TADs.
    """

    def __init__(self, i_file, t_file) -> None:
        if not isinstance(i_file, str) or not os.path.isfile(i_file):
            raise ValueError("Invalid file {}".format(i_file))

        self._figure_size = 10  # What does that mean?

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
                    # for key, d_inter_list in self._d_inter_by_chrom_dict.items():
            # print('\t' + key + ': ' + '{:,}'.format(len(d_inter_list)) + ' interactions')
        print('\t[INFO] ... done.')

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
        # for key, tad_list in self._tads_by_chrom_dict.items():
        # print('\t' + key + ': ' + '{:,}'.format(len(tad_list)) + ' TADs')
        print('\t[INFO] ... done.')
        print('[INFO] ... done.')

    def extract_interactions(self, chrom, begin, end):
        """
        Create reduced interaction list that only contains interactions within the region to be visualized.
        """

        # Check arguments
        if not isinstance(chrom, str):
            raise ValueError("chrom must be a string")
        if not isinstance(begin, int) or not isinstance(end, int):
            raise ValueError("begin and end must be integers")

        # If there are no interactions for this chromosome, return empty list
        if chrom not in self._d_inter_by_chrom_dict:
            return inter_list

        # Extract interactions
        inter_list = []
        for d_inter in self._d_inter_by_chrom_dict[chrom]:

            # Only keep interactions completely within the region we want to visualize
            if d_inter.chrA == chrom and begin < d_inter.fromA and d_inter.toB < end:
                inter_list.append(d_inter)

        print("[INFO] Extracted {} interactions in range {}:{}-{}".format(len(inter_list), chrom, begin, end))
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

    def interaction_to_polygon(self, d_inter, pp_color, pp_alpha):
        """
        Creates a PolygonPatch for a given interaction.
        """
        if not isinstance(d_inter, DiachromaticInteraction):
            raise ValueError("Not a DiachromaticInteraction")

        TANGENT_45 = math.tan(math.pi / 4)  # The lines go up at 45 degree angles
        # a,b,c,d are coordinates of the two digests on the X axis
        a = self.pos_to_coordinate(d_inter.fromA)
        b = self.pos_to_coordinate(d_inter.toA)
        c = self.pos_to_coordinate(d_inter.fromB)
        d = self.pos_to_coordinate(d_inter.toB)
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

        return PolygonPatch(polygon=poly, color=pp_color, alpha=pp_alpha)

    def create_visualization(self,
                             chrom: str,
                             begin: int,
                             end: int,
                             inter_cat_list: List = ['DIX', 'DI', 'UIR', 'UI'],
                             enr_cat_list: List = ['NE', 'EN'],
                             ht_tag_list: List = ['01', '02', '03', '12', '13', '23'],
                             color_i_cats=True):
        """
        Create a TADviz plot.
        """

        i_cat_colors = {
            'DIX': 'orangered',
            'DI': 'orange',
            'UIR': 'green',
            'UI': 'gray'
        }

        ht_tag_colors = {
            '01': 'pink',
            '02': 'red',
            '03': 'lime',
            '12': 'magenta',
            '13': 'blue',
            '23': 'turquoise'
        }

        figure_height = 8
        figure_width = 2 * figure_height

        # Get list of all interactions completely within the region to be visualized
        inter_list = self.extract_interactions(chrom=chrom,
                                               begin=begin,
                                               end=end)

        # If there are no interaction within this region, do nothing
        if len(inter_list) == 0:
            return

        # Get TAD regions for this chromosome
        tads = self._tads_by_chrom_dict[chrom]

        # Get list of interactions to be visualized and collect read pair counts
        d_inter_list = []
        rp_total_list = []
        max_i_dist = 0
        for d_inter in inter_list:
            i_cat = d_inter.get_category()
            e_cat = d_inter.enrichment_status_tag_pair
            ht_tag = d_inter.get_ht_tag()
            # if d_inter.get_category() in inter_cat_list: # Filter for interaction category
            if i_cat in inter_cat_list and e_cat in enr_cat_list and ht_tag in ht_tag_list:
                d_inter_list.append(d_inter)
                if max_i_dist < d_inter.i_dist:
                    max_i_dist = d_inter.i_dist
                rp_total_list.append(d_inter.rp_total)  # (for all interactions or visualized interactions only?)
                # rp_total_list.append(d_inter._log10_pval)

        # Determine read pair counts for quantiles 0.1, ..., 1.0
        quantile_range = np.arange(0.1, 1.1, 0.1)
        quantile_values = np.quantile(rp_total_list, quantile_range)
        print('[INFO] Read pair count quantiles:')
        print('\t[INFO] ' + str(str(quantile_range)))
        print('\t[INFO] ' + str(quantile_values))

        # Get scaling factor
        self._min = begin
        self._max = end
        self._span = end - begin
        self._factor = self._figure_size / self._span

        print("[INFO] Got {} interactions".format(len(d_inter_list)))
        print('[INFO] begin: ' + str(begin))
        print('[INFO] end: ' + str(end))
        print('[TEST] pos_to_coordinate(begin): ' + str(self.pos_to_coordinate(begin)))
        print('[TEST] pos_to_coordinate(end): ' + str(self.pos_to_coordinate(end)))
        print('[TEST] coordinate_to_pos(pos_to_coordinate(begin)): ' + str(
            self.coordinate_to_pos(self.pos_to_coordinate(begin))))
        print('[TEST] coordinate_to_pos(pos_to_coordinate(end)): ' + str(
            self.coordinate_to_pos(self.pos_to_coordinate(end))))

        # Created figure
        fig = plt.figure(1, figsize=(figure_width, figure_height))
        ax = fig.add_subplot(111)

        ax.set_title(str(inter_cat_list) + str(enr_cat_list) + str(ht_tag_list), loc='left', fontsize='x-large')
        ax.set_xlabel('Genomic coordinate', labelpad=12, fontsize='x-large')
        ax.set_ylabel('Interaction distance', labelpad=12, fontsize='x-large')

        # Set limits of x and y axes
        xrange = [self.pos_to_coordinate(begin), self.pos_to_coordinate(end)]
        yrange = [0, 0.5 * (self.pos_to_coordinate(end) - self.pos_to_coordinate(begin))]
        # Uncomment the following line in order to set the upper ylim to the maximum interaction distance
        # yrange = [0, 0.5 * (self.pos_to_coordinate(2*max_i_dist) - self.pos_to_coordinate(max_i_dist))]
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)

        ## Plot the TADs
        for (left_boundary, right_boundary) in tads:
            polypatch = self.tad_to_grey_triangle(start=left_boundary, end=right_boundary)
            ax.add_patch(polypatch)

        # Plot the interaction polygons
        for d_inter in d_inter_list:
            i_cat = d_inter.get_category()
            e_cat = d_inter.enrichment_status_tag_pair
            ht_tag = d_inter.get_ht_tag()
            if i_cat in inter_cat_list and e_cat in enr_cat_list and ht_tag in ht_tag_list:

                # Determine transparency depending on read pair count
                pp_alpha = 0.1
                for i in range(0, 10):
                    if quantile_values[i] <= d_inter.rp_total:
                        pp_alpha = quantile_range[i]
                # if pp_alpha >= 0.8:
                # print(str(d_inter.rp_total) + '\t' + str(pp_alpha))

                # Choose color depending on the interaction category or ht_tag
                if color_i_cats:
                    pp_color = i_cat_colors[i_cat]
                else:
                    pp_color = ht_tag_colors[d_inter.get_ht_tag()]

                # Get and plot interaction polypatch
                polypatch = self.interaction_to_polygon(d_inter=d_inter,
                                                        pp_color=pp_color,
                                                        pp_alpha=pp_alpha)
                ax.add_patch(polypatch)

        # Plot black triangles
        polypatch = self.black_triangle_left(start=begin, end=end)
        ax.add_patch(polypatch)
        polypatch = self.black_triangle_right(start=begin, end=end)
        ax.add_patch(polypatch)

        # Add genomic coordinate labels to x-axis
        ax.set_xticks(ax.get_xticks())
        xtick_labels = ['{:,}'.format(round(self.coordinate_to_pos(x))) for x in ax.get_xticks()]
        ax.set_xticklabels(xtick_labels)

        # Add interaction distance labels to y-axis (correct?)
        ax.set_yticks(ax.get_yticks())
        ytick_labels = ['{:,}'.format(2 * round(self.coordinate_to_pos(y) - begin)) for y in ax.get_yticks()]
        ytick_labels[0] = ''
        ax.set_yticklabels(ytick_labels)

        ax.set_aspect(1)
        plt.show()

