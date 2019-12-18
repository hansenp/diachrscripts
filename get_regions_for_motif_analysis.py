"""
This script is to prepare BED files containing interacting digests and promoter regions used for motif analysis.
It iterates a file with interactions, gene symbols,
strand and TSS information, etc. created with the script 'get_gene_symbols_of_interactions_script.py'.

Six BED files are created:

   1. BED files with regions of interacting digests

      a. of directed interactions ('S' or 'T' in column 3 of file with gene symbols and interactions)
         '<OUT_PREFIX>_directed_digests.bed'

      b. of undirected reference interactions ('URAA' in column 3 of file with gene symbols and interactions)
         '<OUT_PREFIX>_undirected_digests.bed'

   2. BED files for promoter regions

      a. on directed interacting digests
         i. and on the '-' strand
            '<OUT_PREFIX>_directed_minus_tss.bed'

         ii. and on the '+' strand
            '<OUT_PREFIX>_directed_plus_tss.bed'

      b. on undirected interacting digests
         i. and on the '-' strand
            '<OUT_PREFIX>_undirected_minus_tss.bed'

         ii. and on the '+' strand
            '<OUT_PREFIX>_undirected_plus_tss.bed'

Use copy these lines to a text file:

chr9:128781853-128788644;chr9:128855420-128861085	66776	S	TBC1D13;KYAT1	12:0	AA	8.32	+/-	chr9:128787253:+,chr9:128787207:+;chr9:128860532:-
chr17:39397580-39405334;chr17:39461300-39463357	55966	T	FBXL20;CDK12	33:11	AA	7.37	-/+	chr17:39402556:-,chr17:39401626:-;chr17:39461486:+
chr3:101509230-101513979;chr3:101842446-101853540	328467	URAA	SENP7;NFKBIZ	6:3	AA	1.37	-/+	chr3:101513212:-,chr3:101513241:-;chr3:101849514:+
chr21:30736306-30747469;chr21:30819004-30819702	71535	URAI	KRTAP21-2;	8:4	AI	1.64	-/-1	chr21:30747259:-;
chr5:74778161-74782279;chr5:74867764-74868683	85485	URII	;	17:6	II	4.05	-1/-1	;
chr5:88676106-88682018;chr5:88927332-88929569	245314	NA	LINC00461,MEF2C-AS2;	5:1	II	2.21	d/-1	chr5:88676298:-,chr5:88678448:-,chr5:88676218:+;

gzip the test file and use this file as input for this script.
Digest and promoter will be extracted only for the first two lines (directed interactions) and the third line (undirected reference interaction).

The BED files with digest regions are used for motif discovery using DREME.
The BED files with promoter regions are used the analysis of motif occurrences relative to TSS using Centrimo.

We use 'bedtools getfasta' for sequence extraction and gsub for repeat masking.
"""

import argparse
import gzip


### Parse command line
######################

parser = argparse.ArgumentParser(description='Extract regions for motif analysis from file with interactions and gene symbols.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-gs-file', help='Interaction file with gene symbols created with \'get_gene_symbols_of_interactions_script.py\'.')
parser.add_argument('--chrom-info-file', help='File with chromosome sizes <CHROMOSOME NAME>\\t<CHROMOSOME SIZE>.')
parser.add_argument('--up-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--down-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--allowed-enrichment-pair-tags', help='Set of allowed enrichment pair tags for digests (\'AA\', \'AI\',\'IA\',\'II\', ...) separated by \';\'.', default="AA")
parser.add_argument('--allowed-strand-pair-tags', help='Set of allowed strand pair tags for digests (\'-/-\', \'+/-\',\'-/+\',\'-/d\', ...) separated by \';\'.', default="-/-;+/+")

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_gs_file = args.interaction_gs_file # central file that contains a lot of information about interactions
chrom_info_file = args.chrom_info_file # needed to avoid promoter coordinates outside chromosome regions
up_dist = int(args.up_dist) # used to determine start coordinates of promoter regions
down_dist = int(args.down_dist) # used to determine end coordinates of promoter regions
allowed_enrichment_pair_tags = args.allowed_enrichment_pair_tags # only interactions with these enrichment pair tags will be used
allowed_strand_pair_tags = args.allowed_strand_pair_tags # only interactions with these strand pair tags will be used

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] --interaction-gs-file: " + interaction_gs_file)
print("\t[INFO] --chrom-info-file: " + chrom_info_file) #
print("\t[INFO] --up-dist: " + str(up_dist))
print("\t[INFO] --down-dist: " + str(down_dist))
print("\t[INFO] --enrichment-strand-pair-tags: " + allowed_enrichment_pair_tags)
print("\t[INFO] --allowed-strand-pair-tags: " + allowed_strand_pair_tags)


### Prepare output files
########################

# convert aallowed pair tag strings into lists
allowed_strand_pair_tags = str(allowed_strand_pair_tags).split(";")
#print(allowed_strand_pair_tags)
allowed_enrichment_pair_tags = str(allowed_enrichment_pair_tags).split(";")
#print(allowed_enrichment_pair_tags)

# create two BED files for regions of directed interactions
file_name_directed_digests = out_prefix + "_directed_digests.bed"
directed_digests_output_bed = open(file_name_directed_digests, 'wt')
file_name_undirected_digests = out_prefix + "_undirected_digests.bed"
undirected_digests_output_bed = open(file_name_undirected_digests, 'wt')

# create two BED files promoters on digests of directed iteractions, one for TSS on the plus and another one for TSS on the minus strand
file_name_directed_minus_tss = out_prefix + "_directed_minus_tss.bed"
directed_minus_tss_output_bed = open(file_name_directed_minus_tss, 'wt')
file_name_directed_plus_tss = out_prefix + "_directed_plus_tss.bed"
directed_plus_tss_output_bed = open(file_name_directed_plus_tss, 'wt')

# create two BED files promoters on digests of undirected iteractions, one for TSS on the plus and another one for TSS on the minus strand
file_name_undirected_minus_tss = out_prefix + "_undirected_minus_tss.bed"
undirected_minus_tss_output_bed = open(file_name_undirected_minus_tss, 'wt')
file_name_undirected_plus_tss = out_prefix + "_undirected_plus_tss.bed"
undirected_plus_tss_output_bed = open(file_name_undirected_plus_tss, 'wt')


### Prepare variables and data structures
#########################################

chrom_sizes = {}
with open(chrom_info_file, 'rt') as fp:
    line = fp.readline()
    while line:
        chr_name = line.split("\t")[0]
        chr_size = line.split("\t")[1]
        chrom_sizes[chr_name] = int(chr_size)
        line = fp.readline()

directed_minus_tss = set()
directed_plus_tss = set()
undirected_minus_tss = set()
undirected_plus_tss = set()


### Iterate interaction file with gene symbols
##############################################
cnt=0
print("[INFO] Iterating interaction file with gene symbols " + interaction_gs_file + " ...")
with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0 # counter to track progress
    line = fp.readline()

    while line:

        line = line.rstrip()

        n_interaction_total += 1
        if n_interaction_total % 100000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        field = line.split("\t")

        if field[5] not in allowed_enrichment_pair_tags: # use only interactions with specified digest enrichment pair tags
            line = fp.readline()
            continue

        if field[7] not in allowed_strand_pair_tags: # use only interactions with specified digest strand pair tags
            line = fp.readline()
            continue

        # get coordinates of interacting digests
        coords = field[0].split(";")
        coord_a = coords[0]
        coord_b = coords[1]
        chr_a = coord_a.split(":")[0]
        chr_b = coord_b.split(":")[0]
        sta_a = int(coord_a.split(":")[1].split("-")[0])
        sta_b = int(coord_b.split(":")[1].split("-")[0])
        end_a = int(coord_a.split(":")[1].split("-")[1])
        end_b = int(coord_b.split(":")[1].split("-")[1])

        # get category of interacting digest
        interaction_category = field[2]

        if interaction_category == "S" or interaction_category == "T":
            cnt += 1
            # get TSS associated with interacting digests
            tss_d1 = field[8].split(";")[0].split(",")
            for tss in tss_d1:
                chromosome = tss.split(":")[0]
                sta = int(tss.split(":")[1]) - up_dist
                strand = tss.split(":")[2]
                if sta < 0:
                    continue
                end = int(tss.split(":")[1]) + down_dist
                if chrom_sizes[chromosome] < end:
                    continue

                if strand == "-":
                    directed_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                elif strand == "+":
                    directed_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!")

            tss_d2 = field[8].split(";")[1].split(",")
            for tss in tss_d2:
                chromosome = tss.split(":")[0]
                sta = int(tss.split(":")[1]) - up_dist
                strand = tss.split(":")[2]
                if sta < 0:
                    continue
                end = int(tss.split(":")[1]) + down_dist
                if chrom_sizes[chromosome] < end:
                    continue

                if strand == "-":
                    directed_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                elif strand == "+":
                    directed_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!")

            # print digest regions to BED file for directed interactions
            directed_digests_output_bed.write(chr_a + "\t" + str(sta_a) + "\t" + str(end_a) + "\t" + str(cnt) + "\n")
            cnt += 1
            directed_digests_output_bed.write(chr_b + "\t" + str(sta_b) + "\t" + str(end_b) + "\t" + str(cnt) + "\n")

        elif interaction_category == "URAA":
            cnt += 1
            # get TSS associated with the first interacting digest
            tss_d1 = field[8].split(";")[0].split(",")
            for tss in tss_d1:
                chromosome = tss.split(":")[0]
                sta = int(tss.split(":")[1]) - up_dist
                strand = tss.split(":")[2]
                if sta < 0:
                    continue
                end = int(tss.split(":")[1]) + down_dist
                if chrom_sizes[chromosome] < end:
                    continue

                if strand == "-":
                    undirected_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                elif strand == "+":
                    undirected_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!")

            # get TSS associated with the second interacting digest
            tss_d2 = field[8].split(";")[1].split(",")
            for tss in tss_d2:
                chromosome = tss.split(":")[0]
                sta = int(tss.split(":")[1]) - up_dist
                strand = tss.split(":")[2]
                if sta < 0:
                    continue
                end = int(tss.split(":")[1]) + down_dist
                if chrom_sizes[chromosome] < end:
                    continue

                if strand == "-":
                    undirected_minus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                elif strand == "+":
                    undirected_plus_tss.add(chromosome + "\t" + str(sta) + "\t" + str(end) + "\t" + strand)
                else:
                    print("Warning: Strand symbol of TSS was neither \'-\' not \'+\'!")

            # write digest regions to BED file for undirected interactions
            undirected_digests_output_bed.write(chr_a + "\t" + str(sta_a) + "\t" + str(end_a) + "\t" + str(cnt) + "\n")
            cnt += 1
            undirected_digests_output_bed.write(chr_b + "\t" + str(sta_b) + "\t" + str(end_b) + "\t" + str(cnt) + "\n")

        line = fp.readline()


### Write promoter sets to file
###############################

#print("Directed/-")
for promoter in directed_minus_tss:
    directed_minus_tss_output_bed.write(promoter + "\n")
    #print(promoter)
#print("Directed/+")
for promoter in directed_plus_tss:
    directed_plus_tss_output_bed.write(promoter + "\n")
    #print(promoter)
#print("Undirected/-")
for promoter in undirected_minus_tss:
    undirected_minus_tss_output_bed.write(promoter + "\n")
    #print(promoter)
#print("Undirected/+")
for promoter in undirected_plus_tss:
    undirected_plus_tss_output_bed.write(promoter + "\n")
    #print(promoter)


### Close output BED files
##########################

directed_digests_output_bed.close()
undirected_digests_output_bed.close()

directed_minus_tss_output_bed.close()
directed_plus_tss_output_bed.close()

undirected_minus_tss_output_bed.close()
undirected_plus_tss_output_bed.close()
