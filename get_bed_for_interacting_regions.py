"""
This script takes as input an enhanced interaction file and determines two digest sets:

   1. Digests that are involved in at least one directed interaction.
   2. Digest that are involved in undirected reference interactions and in no directed interaction.

The enhanced interaction file is then iterated a second time and the digest regions of each interaction are written
to one of two BED files that:

   1. Contain digest regions that are involved in at least one directed interaction.
   2. Contain digest regions that are involved in undirected reference interactions and in no directed interaction.

Since one digest can be involved in more than one interaction, certain regions may occur in the BED file multiple times.
"""


########################################
# Get BED file with interacting digests

bed_stream_dir_inter_digests_output = open(out_prefix + "_dir_inter_digests.bed", 'wt')
bed_stream_undir_inter_digests_output = open(out_prefix + "_undir_inter_digests.bed", 'wt')

print("[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        x = line.split('\t')
        n_simple = x[4].split(':')[0]
        n_twisted = x[4].split(':')[1]

        if interaction_category == "DIAA" or interaction_category == "DIAI" or interaction_category == "DIII":

            if n_simple < n_twisted:
                bed_name = "T"
            else:
                bed_name = "S"

            # Add digest to set for directed interactions
            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) in digests_from_dir_inter:
                bed_stream_dir_inter_digests_output.write(chr_a + '\t' + str(sta_a) + '\t' + str(end_a) + '\t' + bed_name + '\n')

            # Add digest to set for undirected interactions
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) in digests_from_dir_inter:
                bed_stream_dir_inter_digests_output.write(chr_b + '\t' + str(sta_b) + '\t' + str(end_b) + '\t' + bed_name + '\n')

        if interaction_category == "UIRAA" or interaction_category == "UIRAI" or interaction_category == "UIRII":

            bed_name = "U"

            # Add digest to set for directed interactions
            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) in digests_from_undir_inter_wo_digests_from_dir_inter:
                bed_stream_undir_inter_digests_output.write(chr_a + '\t' + str(sta_a) + '\t' + str(end_a) + '\t' + bed_name + '\n')

            # Add digest to set for undirected interactions
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) in digests_from_undir_inter_wo_digests_from_dir_inter:
                bed_stream_undir_inter_digests_output.write(chr_b + '\t' + str(sta_b) + '\t' + str(end_b) + '\t' + bed_name + '\n')

        line = fp.readline()

bed_stream_dir_inter_digests_output.close()
bed_stream_undir_inter_digests_output.close()