run_for_one_cell_type() {

CELL_TYPE=$1
IN_DIR="results/06_select_reference_interactions/${CELL_TYPE}_RLR"
IN_PREFIX="${CELL_TYPE}_RLR"
OUT_DIR="results/07_analyze_interaction_distances/${CELL_TYPE}_RLR"
OUT_PREFIX="${CELL_TYPE}\_RLR"

mkdir -p $OUT_DIR

# Use Python script in order to get distance vectors for different categories
./05_analyze_interaction_distances.py \
--out-prefix \
$OUT_DIR/$OUT_PREFIX \
--enhanced-interaction-file \
$IN_DIR/${IN_PREFIX}_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz

## Use R-script in order to compare distances of DI, UIR and UI interactions
Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances.r \
$OUT_DIR/ \
$CELL_TYPE\_RLR \
$OUT_DIR/${OUT_PREFIX}_di_ee_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ne_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_en_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_nn_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ee_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ne_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_en_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_nn_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ee_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ne_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_en_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_nn_dist_array.tab

## Use another R-script in order to compare distances of simple and twisted interactions
Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances_st.r \
$OUT_DIR/ \
$CELL_TYPE\_RLR \
$OUT_DIR/${OUT_PREFIX}_di_ee_s_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ne_s_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_en_s_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_nn_s_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ee_t_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ne_t_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_en_t_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_nn_t_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ee_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ne_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_en_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_nn_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ee_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_ne_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_en_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_di_nn_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ee_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ne_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_en_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_nn_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ee_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_ne_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_en_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_uir_nn_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ee_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ne_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_en_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_nn_s_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ee_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_ne_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_en_t_rp_dist_array.tab \
$OUT_DIR/${OUT_PREFIX}_ui_nn_t_rp_dist_array.tab

rm $OUT_DIR/*tab

}

mkdir -p results/07_analyze_interaction_distances/

# Run analysis for individual cell types
#run_for_one_cell_type "MK"
#run_for_one_cell_type "ERY"
#run_for_one_cell_type "NEU"
#run_for_one_cell_type "MON"
#run_for_one_cell_type "MAC_M0"
#run_for_one_cell_type "MAC_M1"
#run_for_one_cell_type "MAC_M2"
#run_for_one_cell_type "EP"
#run_for_one_cell_type "NB"
#run_for_one_cell_type "TB"
#run_for_one_cell_type "FOET"
#run_for_one_cell_type "NCD4"
#run_for_one_cell_type "TCD4"
#run_for_one_cell_type "NACD4"
#run_for_one_cell_type "ACD4"
#run_for_one_cell_type "NCD8"
#run_for_one_cell_type "TCD8"

# Summary - Compare directed and undirected interactions at the level of interactions
#Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats.R \

# Summary - Compare directed and undirected interactions at the level of interactions
#Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats_st.R \

