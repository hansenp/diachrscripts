from diachr import TadBoundarySet, DiachromaticInteraction

mk_tad_boundaries = 'additional_files/javierre_2016/tad_regions_hg38/hglft_genome_TADs_MK_hg38.bed'

tbs = TadBoundarySet(mk_tad_boundaries)

nearest_tad_pos = tbs.get_nearest_tad_boundary('chr1', 10600000)
print(nearest_tad_pos)
nearest_tad_dist = tbs.get_distance_to_nearest_tad_boundary('chr1', 10600000)
print(nearest_tad_dist)


# ## get number of overlapped TAD boundaries for a Digest Pair
# ## chr1    4664940 5492440
# ##chr1    6692440 6942440
# ##chr1    7662440 7757440
# ##chr1    7762440 7958690
#
# ## Assume we have an interaction that spans these digests
# #digest1 = chr1:11645272-11657290
# #digest2 = chr1:11804847-11812432
#
# ia = DiachromaticInteraction(chrA='chr1', fromA=11645272, toA=11657290, statusA='E', chrB='chr1', fromB=11804847, toB=11812432, statusB='N',simple_1=10, simple_2=10, twisted_1=10, twisted_2=10)
#
# tad_boundary_count = tbs.get_overlapped_boundary_count_from_interaction(ia)
#
# print("Tad bounday {}".format(tad_boundary_count))

