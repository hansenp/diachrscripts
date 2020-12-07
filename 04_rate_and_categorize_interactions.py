"""
In this script, the P-values for directionality of interactions are calculated and, based on these P-values,
the interactions are categorized into directed (DI) und undirected interactions (UI).
Undirected reference interactions (UIR) are then selected from UI, which are comparable to DI with respect to
to the total number of read pairs (n).

The input consists of a file in Diachromatic interaction format and a P-value threshold that was determined using our
FDR procedure. The output is again a file in Diachromatic interaction format, but there are two additional columns
on the right for the calculated P-values and interaction categories.
"""

from diachr.binomial_model import BinomialModel

# Parse command line
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('--p-value-threshold', help='P-value threshold for directed interactions.', default=0.001)

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_files_path = args.interaction_files_path
required_replicates = int(args.required_replicates)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --p-value-threshold: " + str(p_value_threshold))


# Determine the time savings acieved by using the dictionary
#p_vaule_factory = BinomialModel()
#p_vaule_factory.measure_time_savings_due_to_pval_dict()

# Load interactions with 'DiachromaticInteractionParser'


# FIRST PASS: CALCULATE P-VALUES AND DEFINE DIRECTED AND UNDIRECTED INTERACTIONS

# SECOND PASS: SELECT UNDIRECTED REFERENCE INTERACTIONS

# WRITE FILE IN DIACHROMATIC INTERACTION FORMAT WITH TWO ADDITIONAL COLUMNS


