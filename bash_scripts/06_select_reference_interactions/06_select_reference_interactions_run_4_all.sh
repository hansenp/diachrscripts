# Select reference interactions for all 17 cell types

mkdir -p results/06_select_reference_interactions/

# MK
mkdir -p results/06_select_reference_interactions/MK
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/MK/MK \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_MK_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# ERY
mkdir -p results/06_select_reference_interactions/ERY
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/ERY/ERY \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_ERY_RALT_0.0018_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# NEU
mkdir -p results/06_select_reference_interactions/NEU
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/NEU/NEU \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_NEU_RALT_0.0011_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# MON
mkdir -p results/06_select_reference_interactions/MON
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/MON/MON \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_MON_RALT_0.0012_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# MAC M0
mkdir -p results/06_select_reference_interactions/MAC_M0
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/MAC_M0/MAC_M0 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_MAC_M0_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# MAC M1
mkdir -p results/06_select_reference_interactions/MAC_M1
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/MAC_M1/MAC_M1 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_MAC_M1_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# MAC M2
mkdir -p results/06_select_reference_interactions/MAC_M2
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/MAC_M2/MAC_M2 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_MAC_M2_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# EP
mkdir -p results/06_select_reference_interactions/EP
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/EP/EP \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_EP_RALT_0.0017_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# NB
mkdir -p results/06_select_reference_interactions/NB
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/NB/NB \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_NB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# TB
mkdir -p results/06_select_reference_interactions/TB
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/TB/TB \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_TB_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# FOET
mkdir -p results/06_select_reference_interactions/FOET
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/FOET/FOET \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_FOET_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# NCD4
mkdir -p results/06_select_reference_interactions/NCD4
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/NCD4/NCD4 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_NCD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# TCD4
mkdir -p results/06_select_reference_interactions/TCD4
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/TCD4/TCD4 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_TCD4_RALT_0.002_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# NACD4
mkdir -p results/06_select_reference_interactions/NACD4
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/NACD4/NACD4 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_NACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# ACD4
mkdir -p results/06_select_reference_interactions/ACD4
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/ACD4/ACD4 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_ACD4_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# NCD8
mkdir -p results/06_select_reference_interactions/NCD8
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/NCD8/NCD8 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_NCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed

# TCD8
mkdir -p results/06_select_reference_interactions/TCD8
./06_select_uir_from_uie_and_uii.py \
--out-prefix \
results/06_select_reference_interactions/TCD8/TCD8 \
--enhanced-interaction-file \
results/05_define_directed_interactions/JAV_TCD8_RALT_0.0019_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz \
--enriched-digests-file \
additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed
