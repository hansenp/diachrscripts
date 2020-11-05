# Perform distance analysis for all 17 cell types
# -----------------------------------------------

# Get path to Python 3.7 binary file
PYTHON_PATH=$1

# Respect left and right, i.e. treat NE and EN as separate categories
#RLR="TRUE"
RLR=$2

# Define input and create directory for output
if [ $RLR = "TRUE" ]; then
  OUT_DIR="results/07_analyze_interaction_distances/RLR"
  declare -a CT_ARRAY=(\
"MK_RLR:results/06_select_reference_interactions/RLR/MK_RLR/MK_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"ERY_RLR:results/06_select_reference_interactions/RLR/ERY_RLR/ERY_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NEU_RLR:results/06_select_reference_interactions/RLR/NEU_RLR/NEU_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MON_RLR:results/06_select_reference_interactions/RLR/MON_RLR/MON_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M0_RLR:results/06_select_reference_interactions/RLR/MAC_M0_RLR/MAC_M0_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M1_RLR:results/06_select_reference_interactions/RLR/MAC_M1_RLR/MAC_M1_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M2_RLR:results/06_select_reference_interactions/RLR/MAC_M2_RLR/MAC_M2_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"EP_RLR:results/06_select_reference_interactions/RLR/EP_RLR/EP_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NB_RLR:results/06_select_reference_interactions/RLR/NB_RLR/NB_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TB_RLR:results/06_select_reference_interactions/RLR/TB_RLR/TB_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"FOET_RLR:results/06_select_reference_interactions/RLR/FOET_RLR/FOET_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NCD4_RLR:results/06_select_reference_interactions/RLR/NCD4_RLR/NCD4_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TCD4_RLR:results/06_select_reference_interactions/RLR/TCD4_RLR/TCD4_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NACD4_RLR:results/06_select_reference_interactions/RLR/NACD4_RLR/NACD4_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"ACD4_RLR:results/06_select_reference_interactions/RLR/ACD4_RLR/ACD4_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NCD8_RLR:results/06_select_reference_interactions/RLR/NCD8_RLR/NCD8_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TCD8_RLR:results/06_select_reference_interactions/RLR/TCD8_RLR/TCD8_RLR_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz")
else
  OUT_DIR="results/07_analyze_interaction_distances"
  declare -a CT_ARRAY=(\
"MK:results/06_select_reference_interactions/MK/MK_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"ERY:results/06_select_reference_interactions/ERY/ERY_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NEU:results/06_select_reference_interactions/NEU/NEU_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MON:results/06_select_reference_interactions/MON/MON_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M0:results/06_select_reference_interactions/MAC_M0/MAC_M0_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M1:results/06_select_reference_interactions/MAC_M1/MAC_M1_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"MAC_M2:results/06_select_reference_interactions/MAC_M2/MAC_M2_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"EP:results/06_select_reference_interactions/EP/EP_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NB:results/06_select_reference_interactions/NB/NB_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TB:results/06_select_reference_interactions/TB/TB_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"FOET:results/06_select_reference_interactions/FOET/FOET_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NCD4:results/06_select_reference_interactions/NCD4/NCD4_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TCD4:results/06_select_reference_interactions/TCD4/TCD4_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NACD4:results/06_select_reference_interactions/NACD4/NACD4_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"ACD4:results/06_select_reference_interactions/ACD4/ACD4_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"NCD8:results/06_select_reference_interactions/NCD8/NCD8_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz" \
"TCD8:results/06_select_reference_interactions/TCD8/TCD8_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz")
fi
mkdir -p $OUT_DIR

# Run interaction distance analysis for all cell types
for ct in "${CT_ARRAY[@]}"; do
  CT_SHORT=$(echo $ct | awk -F: '{print $1}')
  CT_IN_FILE=$(echo $ct | awk -F: '{print $2}')
  CT_OUT_PREFIX=$CT_SHORT
  mkdir -p $OUT_DIR/$CT_OUT_PREFIX

  # Use Python script in order to get distance vectors for different categories
  $PYTHON_PATH ./07_analyze_interaction_distances.py --out-prefix $OUT_DIR/$CT_OUT_PREFIX/$CT_OUT_PREFIX --enhanced-interaction-file $CT_IN_FILE

  ## Use R-script in order to compare distances of DI, UIR and UI interactions
  Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances.r \
  $OUT_DIR/$CT_OUT_PREFIX/ \
  $CT_OUT_PREFIX \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ee_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ne_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_en_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_nn_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ee_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ne_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_en_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_nn_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ee_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ne_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_en_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_nn_dist_array.tab

  ## Use another R-script in order to compare distances of simple and twisted interactions
  Rscript --vanilla rscripts/07_analyze_interaction_distances/interaction_distances_st.r \
  $OUT_DIR/$CT_OUT_PREFIX/ \
  $CT_OUT_PREFIX \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ee_s_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ne_s_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_en_s_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_nn_s_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ee_t_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ne_t_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_en_t_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_nn_t_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ee_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ne_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_en_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_nn_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ee_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_ne_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_en_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_di_nn_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ee_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ne_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_en_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_nn_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ee_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_ne_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_en_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_uir_nn_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ee_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ne_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_en_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_nn_s_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ee_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_ne_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_en_t_rp_dist_array.tab \
  $OUT_DIR/$CT_OUT_PREFIX/${CT_OUT_PREFIX}_ui_nn_t_rp_dist_array.tab

  # Remove files with distances because they are large
  rm $OUT_DIR/$CT_OUT_PREFIX/*tab
done

if [ $RLR = "TRUE" ]; then

  # Summary - Compare directed and undirected interactions at the level of interactions
  Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats.R \
  $OUT_DIR/ \
  "RLR_" \
  " (RLR)" \
  "$OUT_DIR/MK_RLR/MK_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ERY_RLR/ERY_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NEU_RLR/NEU_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MON_RLR/MON_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M0_RLR/MAC_M0_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M1_RLR/MAC_M1_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M2_RLR/MAC_M2_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/EP_RLR/EP_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NB_RLR/NB_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TB_RLR/TB_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/FOET_RLR/FOET_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD4_RLR/NCD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD4_RLR/TCD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NACD4_RLR/NACD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ACD4_RLR/ACD4_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD8_RLR/NCD8_RLR_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD8_RLR/TCD8_RLR_i_distance_statistics_ee_ne_en_nn.tsv"

  # Summary - Compare simple and twisted interactions at the level of interactions and read pairs
  Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats_st.R \
  $OUT_DIR/ \
  "RLR_" \
  " (RLR)" \
  "$OUT_DIR/MK_RLR/MK_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ERY_RLR/ERY_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NEU_RLR/NEU_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MON_RLR/MON_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M0_RLR/MAC_M0_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M1_RLR/MAC_M1_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M2_RLR/MAC_M2_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/EP_RLR/EP_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NB_RLR/NB_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TB_RLR/TB_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/FOET_RLR/FOET_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD4_RLR/NCD4_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD4_RLR/TCD4_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NACD4_RLR/NACD4_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ACD4_RLR/ACD4_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD8_RLR/NCD8_RLR_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD8_RLR/TCD8_RLR_st_distance_statistics_ee_ne_en_nn.tsv"

else

  # Summary - Compare directed and undirected interactions at the level of interactions
  Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats.R \
  $OUT_DIR/ \
  "" \
  "" \
  "$OUT_DIR/MK/MK_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ERY/ERY_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NEU/NEU_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MON/MON_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M0/MAC_M0_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M1/MAC_M1_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M2/MAC_M2_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/EP/EP_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NB/NB_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TB/TB_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/FOET/FOET_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD4/NCD4_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD4/TCD4_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NACD4/NACD4_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ACD4/ACD4_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD8/NCD8_i_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD8/TCD8_i_distance_statistics_ee_ne_en_nn.tsv"

  # Summary - Compare simple and twisted interactions at the level of interactions and read pairs
  Rscript --vanilla rscripts/07_analyze_interaction_distances/analyze_summary_stats_st.R \
  $OUT_DIR/ \
  "" \
  "" \
  "$OUT_DIR/MK/MK_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ERY/ERY_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NEU/NEU_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MON/MON_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M0/MAC_M0_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M1/MAC_M1_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/MAC_M2/MAC_M2_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/EP/EP_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NB/NB_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TB/TB_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/FOET/FOET_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD4/NCD4_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD4/TCD4_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NACD4/NACD4_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/ACD4/ACD4_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/NCD8/NCD8_st_distance_statistics_ee_ne_en_nn.tsv" \
  "$OUT_DIR/TCD8/TCD8_st_distance_statistics_ee_ne_en_nn.tsv"

fi
