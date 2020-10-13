source("rscripts/07_analyze_interaction_distances/interaction_distances_lib.r")

OUT_DIR <- "results/07_analyze_interaction_distances/"

# Create table for all 17 cell types
MK_TAB <- read.table("results/07_analyze_interaction_distances/MK/MK_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
ERY_TAB <- read.table("results/07_analyze_interaction_distances/ERY/ERY_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
NEU_TAB <- read.table("results/07_analyze_interaction_distances/NEU/NEU_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
MON_TAB <- read.table("results/07_analyze_interaction_distances/MON/MON_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
MAC_M0_TAB <- read.table("results/07_analyze_interaction_distances/MAC_M0/MAC_M0_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
MAC_M1_TAB <- read.table("results/07_analyze_interaction_distances/MAC_M1/MAC_M1_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
MAC_M2_TAB <- read.table("results/07_analyze_interaction_distances/MAC_M2/MAC_M2_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
EP_TAB <- read.table("results/07_analyze_interaction_distances/EP/EP_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
NB_TAB <- read.table("results/07_analyze_interaction_distances/NB/NB_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
TB_TAB <- read.table("results/07_analyze_interaction_distances/TB/TB_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
FOET_TAB <- read.table("results/07_analyze_interaction_distances/FOET/FOET_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
NCD4_TAB <- read.table("results/07_analyze_interaction_distances/NCD4/NCD4_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
TCD4_TAB <- read.table("results/07_analyze_interaction_distances/TCD4/TCD4_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
NACD4_TAB <- read.table("results/07_analyze_interaction_distances/NACD4/NACD4_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
ACD4_TAB <- read.table("results/07_analyze_interaction_distances/ACD4/ACD4_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
NCD8_TAB <- read.table("results/07_analyze_interaction_distances/NCD8/NCD8_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)
TCD8_TAB <- read.table("results/07_analyze_interaction_distances/TCD8/TCD8_i_distance_statistics_ee_ne_en_nn.tsv", header=TRUE)

ALL_TAB <- rbind(
  MK_TAB,
  ERY_TAB,
  NEU_TAB,
  MON_TAB,
  MAC_M0_TAB,
  MAC_M1_TAB,
  MAC_M2_TAB,
  EP_TAB,
  NB_TAB,
  TB_TAB,
  FOET_TAB,
  NCD4_TAB,
  TCD4_TAB,
  NACD4_TAB,
  ACD4_TAB,
  NCD8_TAB,
  TCD8_TAB
)

cairo_pdf(paste(OUT_DIR,"interaction_distance_summary_stats.pdf",sep=""), width=18, height=56)
par(mfrow=c(14,3), oma = c(0, 0, 2, 0))

# Interaction numbers - Overview
N_TAB_DI <- ALL_TAB[,c("DI_EE_N","DI_NE_N","DI_EN_N","DI_NN_N")]
boxplot(
  N_TAB_DI,
  main="Directed interaction numbers for the 17 cell types",
  col=directed_color,
  names=c("EE", "NE", "EN", "NN")
  )

N_TAB_UIR <- ALL_TAB[,c("UIR_EE_N","UIR_NE_N","UIR_EN_N","UIR_NN_N")]
boxplot(
  N_TAB_UIR,
  main="Undirected reference interaction numbers for the 17 cell types",
  col=undirected_ref_color,
  names=c("EE", "NE", "EN", "NN")
)


N_TAB_UI <- ALL_TAB[,c("UI_EE_N","UI_NE_N","UI_EN_N","UI_NN_N")]
boxplot(
  N_TAB_UI,
  main="Undirected interaction numbers for the 17 cell types",
  col=undirected_color,
  names=c("EE", "NE", "EN", "NN")
)



# Interaction numbers - Differences of NE and EN

N_DIFF_NE_EN_DI <- ALL_TAB[,"DI_NE_N"]-ALL_TAB[,"DI_EN_N"]
N_DIFF_NE_EN_UIR <- ALL_TAB[,"UIR_NE_N"]-ALL_TAB[,"UIR_EN_N"]
N_DIFF_NE_EN_UI <- ALL_TAB[,"UI_NE_N"]-ALL_TAB[,"UI_EN_N"]

YMAX <- max(N_DIFF_NE_EN_DI,N_DIFF_NE_EN_UIR,N_DIFF_NE_EN_UI) + 1500
YMIN <- min(N_DIFF_NE_EN_DI,N_DIFF_NE_EN_UIR,N_DIFF_NE_EN_UI)

N_DIFF_NE_EN_DI_TT <- t.test(N_DIFF_NE_EN_DI)
N_DIFF_NE_EN_UIR_TT <- t.test(N_DIFF_NE_EN_UIR)
N_DIFF_NE_EN_UI_TT <- t.test(N_DIFF_NE_EN_UI)

boxplot(
  cbind(N_DIFF_NE_EN_DI,N_DIFF_NE_EN_UIR,N_DIFF_NE_EN_UI),
  main="Differences of NE and EN interactions",
  col=c(directed_color,undirected_ref_color,undirected_color),
  names=c("DI", "UIR", "UI"),
  ylim=c(YMIN,YMAX)
  )
abline(h=0,col="grey")

v_space <- 1000
text(c(1:3),
     c(max(N_DIFF_NE_EN_DI)+v_space,
       max(N_DIFF_NE_EN_UIR)+v_space,
       max(N_DIFF_NE_EN_UI)+v_space),
     c(paste("p = ", formatC(N_DIFF_NE_EN_DI_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(N_DIFF_NE_EN_UIR_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(N_DIFF_NE_EN_UI_TT$p.value, format = "e", digits = 2), sep="")
        )
      )

plot.new()
plot.new()

# Median interaction distances - Overview
# ---------------------------------------

MED_TAB_DI <- ALL_TAB[,c("DI_EE_MED","DI_NE_MED","DI_EN_MED","DI_NN_MED")]
boxplot(
  MED_TAB_DI,
  main="Median distances - Directed interactions",
  col=directed_color,
  names=c("EE", "NE", "EN", "NN")
)

MED_TAB_UIR <- ALL_TAB[,c("UIR_EE_MED","UIR_NE_MED","UIR_EN_MED","UIR_NN_MED")]
boxplot(
  MED_TAB_UIR,
  main="Median distances - Undirected reference interactions",
  col=undirected_ref_color,
  names=c("EE", "NE", "EN", "NN")
)

MED_TAB_UI <- ALL_TAB[,c("UI_EE_MED","UI_NE_MED","UI_EN_MED","UI_NN_MED")]
boxplot(
  MED_TAB_UI,
  main="Median distances - Undirected interactions",
  col=undirected_color,
  names=c("EE", "NE", "EN", "NN")
)

MED_TAB_DI_UIR_NE_EN <- ALL_TAB[,c("DI_EE_MED","UIR_EE_MED","DI_NE_MED","UIR_NE_MED","DI_NE_MED","UIR_NE_MED","DI_NN_MED","UIR_NN_MED")]
boxplot(
  MED_TAB_DI_UIR_NE_EN,
  main="Median distances - Directed and undirected reference interactions",
  col=c(directed_color,
        undirected_ref_color
        ),
  names=c("DI-EE", "UIR-EE", "DI-NE", "UIR-NE", "DI-EN", "UIR-EN", "DI-NN", "UIR-NN")
)


# Median interaction distances - Differences of NE and EN
# -------------------------------------------------------

MED_DIFF_NE_EN_DI <- ALL_TAB[,"DI_NE_MED"]-ALL_TAB[,"DI_EN_MED"]
MED_DIFF_NE_EN_UIR <- ALL_TAB[,"UIR_NE_MED"]-ALL_TAB[,"UIR_EN_MED"]
MED_DIFF_NE_EN_UI <- ALL_TAB[,"UI_NE_MED"]-ALL_TAB[,"UI_EN_MED"]

YMAX <- max(MED_DIFF_NE_EN_DI,MED_DIFF_NE_EN_UIR,MED_DIFF_NE_EN_UI) + 1500
YMIN <- min(MED_DIFF_NE_EN_DI,MED_DIFF_NE_EN_UIR,MED_DIFF_NE_EN_UI)

MED_DIFF_NE_EN_DI_TT <- t.test(MED_DIFF_NE_EN_DI)
MED_DIFF_NE_EN_UIR_TT <- t.test(MED_DIFF_NE_EN_UIR)
MED_DIFF_NE_EN_UI_TT <- t.test(MED_DIFF_NE_EN_UI)

boxplot(
  cbind(MED_DIFF_NE_EN_DI,MED_DIFF_NE_EN_UIR,MED_DIFF_NE_EN_UI),
  main="Differences of NE and EN median interaction distances",
  col=c(directed_color,undirected_ref_color,undirected_color),
  names=c("DI", "UIR", "UI"),
  ylim=c(YMIN,YMAX)
)
abline(h=0,col="grey")

v_space <- 1000
text(c(1:3),
     c(max(MED_DIFF_NE_EN_DI)+v_space,
       max(MED_DIFF_NE_EN_UIR)+v_space,
       max(MED_DIFF_NE_EN_UI)+v_space),
     c(paste("p = ", formatC(MED_DIFF_NE_EN_DI_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(MED_DIFF_NE_EN_UIR_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(MED_DIFF_NE_EN_UI_TT$p.value, format = "e", digits = 2), sep="")
     )
)


plot.new()

# Median absolute deviations - Overview
# -------------------------------------

IQR_TAB_DI <- ALL_TAB[,c("DI_EE_IQR","DI_NE_IQR","DI_EN_IQR","DI_NN_IQR")]
boxplot(
  IQR_TAB_DI,
  main="IQR distances - Directed interactions",
  col=directed_color,
  names=c("EE", "NE", "EN", "NN")
)

IQR_TAB_UIR <- ALL_TAB[,c("UIR_EE_IQR","UIR_NE_IQR","UIR_EN_IQR","UIR_NN_IQR")]
boxplot(
  IQR_TAB_UIR,
  main="IQR - Undirected reference interactions",
  col=undirected_ref_color,
  names=c("EE", "NE", "EN", "NN")
)

IQR_TAB_UI <- ALL_TAB[,c("UI_EE_IQR","UI_NE_IQR","UI_EN_IQR","UI_NN_IQR")]
boxplot(
  IQR_TAB_UI,
  main="IQR - Undirected interactions",
  col=undirected_color,
  names=c("EE", "NE", "EN", "NN")
)

IQR_TAB_DI_UIR_NE_EN <- ALL_TAB[,c("DI_EE_IQR","UIR_EE_IQR","DI_NE_IQR","UIR_NE_IQR","DI_NE_IQR","UIR_NE_IQR","DI_NN_IQR","UIR_NN_IQR")]
boxplot(
  IQR_TAB_DI_UIR_NE_EN,
  main="IQR distances - Directed and undirected reference interactions",
  col=c(directed_color,
        undirected_ref_color
  ),
  names=c("DI-EE", "UIR-EE", "DI-NE", "UIR-NE", "DI-EN", "UIR-EN", "DI-NN", "UIR-NN")
)


# IQR - Differences of NE and EN
# -------------------------------------------------------

IQR_DIFF_NE_EN_DI <- ALL_TAB[,"DI_NE_IQR"]-ALL_TAB[,"DI_EN_IQR"]
IQR_DIFF_NE_EN_UIR <- ALL_TAB[,"UIR_NE_IQR"]-ALL_TAB[,"UIR_EN_IQR"]
IQR_DIFF_NE_EN_UI <- ALL_TAB[,"UI_NE_IQR"]-ALL_TAB[,"UI_EN_IQR"]

YMAX <- max(IQR_DIFF_NE_EN_DI,IQR_DIFF_NE_EN_UIR,IQR_DIFF_NE_EN_UI) + 1500
YMIN <- min(IQR_DIFF_NE_EN_DI,IQR_DIFF_NE_EN_UIR,IQR_DIFF_NE_EN_UI)

IQR_DIFF_NE_EN_DI_TT <- t.test(IQR_DIFF_NE_EN_DI)
IQR_DIFF_NE_EN_UIR_TT <- t.test(IQR_DIFF_NE_EN_UIR)
IQR_DIFF_NE_EN_UI_TT <- t.test(IQR_DIFF_NE_EN_UI)

boxplot(
  cbind(IQR_DIFF_NE_EN_DI,IQR_DIFF_NE_EN_UIR,IQR_DIFF_NE_EN_UI),
  main="Differences of NE and EN - IQR interaction distances",
  col=c(directed_color,undirected_ref_color,undirected_color),
  names=c("DI", "UIR", "UI"),
  ylim=c(YMIN,YMAX)
)
abline(h=0,col="grey")

v_space <- 1000
text(c(1:3),
     c(max(IQR_DIFF_NE_EN_DI)+v_space,
       max(IQR_DIFF_NE_EN_UIR)+v_space,
       max(IQR_DIFF_NE_EN_UI)+v_space),
     c(paste("p = ", formatC(IQR_DIFF_NE_EN_DI_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(IQR_DIFF_NE_EN_UIR_TT$p.value, format = "e", digits = 2), sep=""),
       paste("p = ", formatC(IQR_DIFF_NE_EN_UI_TT$p.value, format = "e", digits = 2), sep="")
     )
)

dev.off()





