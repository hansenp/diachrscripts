import diachrscripts_toolkit as dclass
import argparse
import gzip


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--up-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--down-dist', help='Number of bases upstream of TSS.', default=1000)

args = parser.parse_args()
out_prefix = args.out_prefix
ref_gene_file = args.ref_gene_file
up_dist = int(args.up_dist)
down_dist = int(args.down_dist)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] refGene file: " + ref_gene_file)
print("\t[INFO] --up-dist: " + str(up_dist))
print("\t[INFO] --down-dist: " + str(down_dist))

# create two BED files, one for the plus and another one for the minus strand
file_name_minus = out_prefix + "_tss_minus.bed"
minus_output_bed = open(file_name_minus, 'wt')
file_name_plus = out_prefix + "_tss_plus.bed"
plus_output_bed = open(file_name_plus, 'wt')

# prepare variables and data structures
ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

# construct BED lines and make TSS unique
coord_set = set()
for coord_key in ref_gene_tss_map.tss_coord_dict.keys():
    if coord_key not in coord_set:
        coord_set.add(coord_key)
        coord_key_arr = coord_key.split(':')
        chromosome = coord_key_arr[0]
        start = int(coord_key_arr[1]) - up_dist
        end = int(coord_key_arr[1]) + down_dist
        strand = ref_gene_tss_map.get_coord_strand(coord_key)
        if 0 < start:
            if strand == '-':
                minus_output_bed.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\n")
            else:
                plus_output_bed.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\n")

plus_output_bed.close()
minus_output_bed.close()