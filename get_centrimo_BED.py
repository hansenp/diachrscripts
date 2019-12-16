import diachrscripts_toolkit as dclass
import argparse
import gzip


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-gs-file', help='Interaction file with gene symbols created with \'get_gene_symbols_of_interactions_script.py\'.')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--up-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--down-dist', help='Number of bases upstream of TSS.', default=1000)

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_gs_file = args.interaction_gs_file
ref_gene_file = args.ref_gene_file
up_dist = int(args.up_dist)
down_dist = int(args.down_dist)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] --interaction-gs-file: " + interaction_gs_file)
print("\t[INFO] --refGene file: " + ref_gene_file)
print("\t[INFO] --up-dist: " + str(up_dist))
print("\t[INFO] --down-dist: " + str(down_dist))


### Prepare output files
########################

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

# create two BED files promoters on digests of directed iteractions, one for TSS on the plus and another one for TSS on the minus strand
file_name_undirected_minus = out_prefix + "_tss_undirected_minus.bed"
undirected_minus_output_bed = open(file_name_undirected_minus, 'wt')
file_name_undirected_plus = out_prefix + "_tss_undirected_plus.bed"
undirected_plus_output_bed = open(file_name_undirected_plus, 'wt')


### Prepare variables and data structures
#########################################

# read refGene.txt.gz to data structure
ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

directed_gene_symbols = set()
undirected_gene_symbols = set()

### Iterate interaction file with gene symbols
##############################################

print("[INFO] Iterating interaction file with gene symbols " + interaction_gs_file + " ...")
with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()

    while line:

        n_interaction_total += 1
        if n_interaction_total % 10000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        # get coordinates of interacting digests
        field = line.split("\t")
        coords = field[0].split(";")
        coord_a = coords[0]
        coord_b = coords[1]
        chr_a = coord_a.split(":")[0]
        chr_b = coord_b.split(":")[0]
        sta_a = coord_a.split(":")[1].split("-")[0]
        sta_b = coord_b.split(":")[1].split("-")[0]
        end_a = coord_a.split(":")[1].split("-")[1]
        end_b = coord_b.split(":")[1].split("-")[1]

        # get category of interacting digest
        interaction_category = field[2]

        if interaction_category == "S" or interaction_category == "T":

            # add gene symbols to set for directed interactions
            symbols_a = field[3].split(";")[0].split(",")
            for sym_a in symbols_a:
                if sym_a != "":
                    directed_gene_symbols.add(sym_a)
            symbols_b = field[3].split(";")[1].split(",")
            for sym_b in symbols_b:
                if sym_b != "":
                    directed_gene_symbols.add(sym_b)

            # print digest regions to BED file for directed interactions
            directed_digests_output_bed.write(chr_a + "\t" + sta_a + "\t" + end_a + "\n")
            directed_digests_output_bed.write(chr_b + "\t" + sta_b + "\t" + end_b + "\n")

        elif interaction_category == "URAA":

            # add gene symbols to set for directed interactions
            symbols_a = field[3].split(";")[0].split(",")
            for sym_a in symbols_a:
                if sym_a != "":
                    undirected_gene_symbols.add(sym_a)
            symbols_b = field[3].split(";")[1].split(",")
            for sym_b in symbols_b:
                if sym_b != "":
                    undirected_gene_symbols.add(sym_b)

            # print digest regions to BED file for directed interactions
            undirected_digests_output_bed.write(chr_a + "\t" + sta_a + "\t" + end_a + "\n")
            undirected_digests_output_bed.write(chr_b + "\t" + sta_b + "\t" + end_b + "\n")

        line = fp.readline()


### Iterate TSS
###############

# construct BED lines and make TSS unique
coord_set = set()
for coord_key in ref_gene_tss_map.tss_coord_dict.keys():

    #  find out if gene symbol(s) is contained in set for directed and undirected digests
    goto_next_coord = True
    for sym in ref_gene_tss_map.get_coord_symbols(coord_key):
        if sym in directed_gene_symbols  or sym in undirected_gene_symbols:
            goto_next_coord = False

    if goto_next_coord:
        print("TSS not in digest.")
        continue
    else:
        print("TSS in digest.")

    if coord_key not in coord_set:
        coord_set.add(coord_key)
        coord_key_arr = coord_key.split(':')
        chromosome = coord_key_arr[0]
        start = int(coord_key_arr[1]) - up_dist
        end = int(coord_key_arr[1]) + down_dist
        strand = ref_gene_tss_map.get_coord_strand(coord_key)
        if 0 < start:
            if strand == '-':
                directed_minus_tss_output_bed.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\n")
            else:
                directed_plus_tss_output_bed.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\n")

directed_digests_output_bed.close()
undirected_digests_output_bed.close()

directed_minus_tss_output_bed.close()
directed_plus_tss_output_bed.close()

undirected_minus_output_bed.close()
undirected_plus_output_bed.close()