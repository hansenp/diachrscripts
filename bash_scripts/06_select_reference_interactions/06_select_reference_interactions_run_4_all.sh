# Select reference interactions for all 17 cell types
# ---------------------------------------------------

# Get path to Python 3.7 binary file
#PYTHON_PATH="/Users/hansep/anaconda2/envs/diachscripts_p37env/bin/python3.7"
PYTHON_PATH=$1

# Respect left and right, i.e. treat NE and EN as separate categories
#RLR="TRUE"
RLR=$2

# Get path to enriched digest file
ENRICHED_DIGEST_FILE="data/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed"

# Create directory for output
if [ $RLR = "TRUE" ]; then
  OUT_DIR="results/06_select_reference_interactions/RLR"
else
  OUT_DIR="results/06_select_reference_interactions_x"
fi
mkdir -p $OUT_DIR

# Abbreviations for cell types and input files
declare -a CT_ARRAY=(
"MK:results/05_define_directed_interactions/JAV_MK_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"ERY:results/05_define_directed_interactions/JAV_ERY_RALT_0.0018_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"NEU:results/05_define_directed_interactions/JAV_NEU_RALT_0.0011_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"MON:results/05_define_directed_interactions/JAV_MON_RALT_0.0012_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"MAC_M0:results/05_define_directed_interactions/JAV_MAC_M0_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"MAC_M1:results/05_define_directed_interactions/JAV_MAC_M1_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"MAC_M2:results/05_define_directed_interactions/JAV_MAC_M2_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"EP:results/05_define_directed_interactions/JAV_EP_RALT_0.0017_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"NB:results/05_define_directed_interactions/JAV_NB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"TB:results/05_define_directed_interactions/JAV_TB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"FOET:results/05_define_directed_interactions/JAV_FOET_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"NCD4:results/05_define_directed_interactions/JAV_NCD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"TCD4:results/05_define_directed_interactions/JAV_TCD4_RALT_0.002_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"NACD4:results/05_define_directed_interactions/JAV_NACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"ACD4:results/05_define_directed_interactions/JAV_ACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"NCD8:results/05_define_directed_interactions/JAV_NCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
"TCD8:results/05_define_directed_interactions/JAV_TCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz")

# Run reference selection for all cell types
for ct in "${CT_ARRAY[@]}"; do
  CT_SHORT=$(echo $ct | awk -F: '{print $1}')
  CT_IN_FILE=$(echo $ct | awk -F: '{print $2}')
  echo $CT_SHORT
  echo $CT_IN_FILE

  if [ $RLR = "TRUE" ]; then
    CT_OUT_PREFIX="${CT_SHORT}_RLR"
    mkdir -p $OUT_DIR/$CT_OUT_PREFIX
    $PYTHON_PATH ./06_select_uir_from_uie_and_uii.py \
    --out-prefix \
    $OUT_DIR/$CT_OUT_PREFIX/$CT_OUT_PREFIX \
    --enhanced-interaction-file \
    $CT_IN_FILE \
    --enriched-digests-file \
    $ENRICHED_DIGEST_FILE \
    --respect-left-right
  else
    CT_OUT_PREFIX=$CT_SHORT
    mkdir -p $OUT_DIR/$CT_OUT_PREFIX
    $PYTHON_PATH ./06_select_uir_from_uie_and_uii.py \
    --out-prefix \
    $OUT_DIR/$CT_OUT_PREFIX/$CT_OUT_PREFIX \
    --enhanced-interaction-file \
    $CT_IN_FILE \
    --enriched-digests-file \
    $ENRICHED_DIGEST_FILE
  fi
done
