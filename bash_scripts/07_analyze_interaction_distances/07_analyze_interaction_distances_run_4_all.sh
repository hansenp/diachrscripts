run_for_one_cell_type() {

CELL_TYPE=$1
IN_DIR="results/06_select_reference_interactions/${CELL_TYPE}_RLR"
IN_PREFIX="${CELL_TYPE}_RLR"
OUT_DIR="results/07_analyze_interaction_distances/RLR/${CELL_TYPE}_RLR"
OUT_PREFIX="${CELL_TYPE}_RLR"

mkdir -p $OUT_DIR

# Use Python script in order to get distance vectors for different categories
./07_analyze_interaction_distances.py \
--out-prefix \
$OUT_DIR/$OUT_PREFIX \
--enhanced-interaction-file \
$IN_DIR/${IN_PREFIX}_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz

## Use R-script in order to compare distances of DI, UIR and UI interactions
Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances.r \
$OUT_DIR/ \
$OUT_PREFIX \
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
OUT_PREFIX \
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
Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats.R \
results/07_analyze_interaction_distances/RLR/ \
"RLR_" \
" (RLR)" \
"results/07_analyze_interaction_distances/RLR/MK_RLR/MK_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/ERY_RLR/ERY_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/NEU_RLR/NEU_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/MON_RLR/MON_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/MAC_M0_RLR/MAC_M0_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/MAC_M1_RLR/MAC_M1_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/MAC_M2_RLR/MAC_M2_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/EP_RLR/EP_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/NB_RLR/NB_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/TB_RLR/TB_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/FOET_RLR/FOET_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/NCD4_RLR/NCD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/TCD4_RLR/TCD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/NACD4_RLR/NACD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/ACD4_RLR/ACD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/NCD8_RLR/NCD8_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/RLR/TCD8_RLR/TCD8_RLR_i_distance_statistics_ee_ne_en_nn.tsv"


Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats.R \
results/07_analyze_interaction_distances/ \
"" \
"" \
"results/07_analyze_interaction_distances/MK/MK_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/ERY/ERY_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/NEU/NEU_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/MON/MON_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/MAC_M0/MAC_M0_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/MAC_M1/MAC_M1_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/MAC_M2/MAC_M2_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/EP/EP_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/NB/NB_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/TB/TB_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/FOET/FOET_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/NCD4/NCD4_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/TCD4/TCD4_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/NACD4/NACD4_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/ACD4/ACD4_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/NCD8/NCD8_i_distance_statistics_ee_ne_en_nn.tsv" \
"results/07_analyze_interaction_distances/TCD8/TCD8_i_distance_statistics_ee_ne_en_nn.tsv"


# Summary - Compare directed and undirected interactions at the level of interactions
#Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats_st.R \

