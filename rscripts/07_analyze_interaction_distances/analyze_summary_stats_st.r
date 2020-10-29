source("rscripts/07_analyze_interaction_distances/interaction_distances_lib.r")

# Read commannd line arguments
args = commandArgs(trailingOnly=TRUE)
for(arg in args){
  print(arg)
}

OUT_DIR <- args[1]
OUT_PREFIX <- args[2]
MM_TITLE_SUFFIX <- args[3]

MK_TAB <- read.table(args[4], header=TRUE)
ERY_TAB <- read.table(args[5], header=TRUE)
NEU_TAB <- read.table(args[6], header=TRUE)
MON_TAB <- read.table(args[7], header=TRUE)
MAC_M0_TAB <- read.table(args[8], header=TRUE)
MAC_M1_TAB <- read.table(args[9], header=TRUE)
MAC_M2_TAB <- read.table(args[10], header=TRUE)
EP_TAB <- read.table(args[11], header=TRUE)
NB_TAB <- read.table(args[12], header=TRUE)
TB_TAB <- read.table(args[13], header=TRUE)
FOET_TAB <- read.table(args[14], header=TRUE)
NCD4_TAB <- read.table(args[15], header=TRUE)
TCD4_TAB <- read.table(args[16], header=TRUE)
NACD4_TAB <- read.table(args[17], header=TRUE)
ACD4_TAB <- read.table(args[18], header=TRUE)
NCD8_TAB <- read.table(args[19], header=TRUE)
TCD8_TAB <- read.table(args[20], header=TRUE)

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

function.get_simple_twisted_boxplot <- function(M_TITLE, YLAB, BOX_COLOR, TAB, YLIMS)
{
  if(YLIMS[1] != FALSE) {
    y_min <- YLIMS[1]
    y_max <- YLIMS[2] 
  }  else {
    y_min <- min(TAB)
    y_max <- max(TAB)   
  }
  
  # Create boxplot for distributions of n
  boxplot(
    TAB,
    main=M_TITLE,
    col=BOX_COLOR,
    xlab="Enrichment status and orientation category",
    names=c("EE-S",
            "EE-T",
            "NE-S",
            "NE-T",
            "EN-S",
            "EN-T",
            "NN-S",
            "NN-T"
    ),
    ylab=YLAB,
    ylim=c(y_min,y_max),
    border=c(simple_color, twisted_color)
  )
}

function.get_simple_twisted_diff_boxplot <- function(M_TITLE, YLAB, BOX_COLOR, TAB, YLIMS)
{
  # Calculate differences between simple and twisted
  ST_DIFF_EE <- TAB[,1]-TAB[,2]
  ST_DIFF_NE <- TAB[,3]-TAB[,4]
  ST_DIFF_EN <- TAB[,5]-TAB[,6]
  ST_DIFF_NN <- TAB[,7]-TAB[,8]
  
  # Perform t-test
  ST_DIFF_EE_TT <- t.test(ST_DIFF_EE)
  ST_DIFF_NE_TT <- t.test(ST_DIFF_NE)
  ST_DIFF_EN_TT <- t.test(ST_DIFF_EN)
  ST_DIFF_NN_TT <- t.test(ST_DIFF_NN)
  
  # Create boxplot for distributions of differences of n between simple and twisted
  if(YLIMS[1] != FALSE) {
    y_min <- YLIMS[1]
    y_max <- YLIMS[2] 
  }  else {
    y_min <- min(c(ST_DIFF_EE,ST_DIFF_NE,ST_DIFF_EN,ST_DIFF_NN))
    y_max <- max(c(ST_DIFF_EE,ST_DIFF_NE,ST_DIFF_EN,ST_DIFF_NN))   
  }
  y_range <- abs(y_min-y_max)
  YMIN <- y_min
  YMAX <- y_max + y_range/10
  
  boxplot(
    cbind(ST_DIFF_EE,ST_DIFF_NE,ST_DIFF_EN,ST_DIFF_NN),
    ylim=c(YMIN,YMAX),
    main=M_TITLE,
    ylab=YLAB,
    col=BOX_COLOR,
    xlab = "Ennrichment status category",
    names=c("EE",
            "NE",
            "EN",
            "NN"
    ),
  )
  abline(h=0,col="grey")
  
  # Add labels with P-values for simple twisted column pairs
  v_space <- y_range/20
  text(c(1:4),
       c(max(ST_DIFF_EE) + v_space,
         max(ST_DIFF_NE) + v_space,
         max(ST_DIFF_EN) + v_space,
         max(ST_DIFF_NN) + v_space),
       c(paste("p = ", formatC(ST_DIFF_EE_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(ST_DIFF_NE_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(ST_DIFF_EN_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(ST_DIFF_NN_TT$p.value, format = "e", digits = 2), sep="")
       )
  )
  
  # Return vector with differennces of simple and twisted
  #return(c(ST_DIFF_EE,ST_DIFF_NE,ST_DIFF_EN,ST_DIFF_NN))
  return(c(ST_DIFF_NE,ST_DIFF_EN))
  
}


# Interaction and read pair numbers
# ---------------------------------
cairo_pdf(paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_n.pdf",sep=""), width=24, height=16)
par(mfrow=c(4,4), oma = c(0, 0, 3, 0))
  
  # Directed interactions
  N_TAB_DI <- ALL_TAB[,c("DI_EE_S_N",
                         "DI_EE_T_N",
                         "DI_NE_S_N",
                         "DI_NE_T_N",
                         "DI_EN_S_N",
                         "DI_EN_T_N",
                         "DI_NN_S_N",
                         "DI_NN_T_N"
  )]
  function.get_simple_twisted_boxplot(
    "Directed simple and twisted interaction numbers",
    "Interaction number",
    directed_color,
    N_TAB_DI,
    FALSE)
  
  # Read pair numbers
  # -----------------
  
  # Directed interactions
  N_TAB_DI_RP <- ALL_TAB[,c("DI_EE_S_RP_N",
                            "DI_EE_T_RP_N",
                            "DI_NE_S_RP_N",
                            "DI_NE_T_RP_N",
                            "DI_EN_S_RP_N",
                            "DI_EN_T_RP_N",
                            "DI_NN_S_RP_N",
                            "DI_NN_T_RP_N"
  )]
  function.get_simple_twisted_boxplot(
    "Simple and twisted read pairs in directed interactions",
    "Interaction number",
    directed_color,
    N_TAB_DI_RP,
    FALSE)
  
  # Undirected reference interactions
  N_TAB_UIR_RP <- ALL_TAB[,c("UIR_EE_S_RP_N",
                             "UIR_EE_T_RP_N",
                             "UIR_NE_S_RP_N",
                             "UIR_NE_T_RP_N",
                             "UIR_EN_S_RP_N",
                             "UIR_EN_T_RP_N",
                             "UIR_NN_S_RP_N",
                             "UIR_NN_T_RP_N"
  )]
  function.get_simple_twisted_boxplot(
    "Simple and twisted read pairs in undirected reference interactions",
    "Interaction number",
    undirected_ref_color,
    N_TAB_UIR_RP,
    FALSE)
  
  # Undirected interactions
  N_TAB_UI_RP <- ALL_TAB[,c("UI_EE_S_RP_N",
                             "UI_EE_T_RP_N",
                             "UI_NE_S_RP_N",
                             "UI_NE_T_RP_N",
                             "UI_EN_S_RP_N",
                             "UI_EN_T_RP_N",
                             "UI_NN_S_RP_N",
                             "UI_NN_T_RP_N"
  )]
  function.get_simple_twisted_boxplot(
    "Simple and twisted read pairs in undirected interactions",
    "Interaction number",
    undirected_color,
    N_TAB_UI_RP,
    FALSE)
  
  YL <- c(
    min(N_TAB_DI_RP[,c("DI_NE_S_RP_N","DI_NE_T_RP_N","DI_EN_S_RP_N","DI_EN_T_RP_N")],
        N_TAB_UIR_RP[,c("UIR_NE_S_RP_N","UIR_NE_T_RP_N","UIR_EN_S_RP_N","UIR_EN_T_RP_N")]),
    max(N_TAB_DI_RP[,c("DI_NE_S_RP_N","DI_NE_T_RP_N","DI_EN_S_RP_N","DI_EN_T_RP_N")],
        N_TAB_UIR_RP[,c("UIR_NE_S_RP_N","UIR_NE_T_RP_N","UIR_EN_S_RP_N","UIR_EN_T_RP_N")])
  )
  
  plot.new()
  
  # Directed interactions
  N_TAB_DI_RP <- ALL_TAB[,c("DI_EE_S_RP_N",
                            "DI_EE_T_RP_N",
                            "DI_NE_S_RP_N",
                            "DI_NE_T_RP_N",
                            "DI_EN_S_RP_N",
                            "DI_EN_T_RP_N",
                            "DI_NN_S_RP_N",
                            "DI_NN_T_RP_N"
  )]
  function.get_simple_twisted_boxplot(
    "Simple and twisted read pairs in directed interactions",
    "Interaction number",
    directed_color,
    N_TAB_DI_RP,
    YL)
  
  # Undirected reference interactions
  N_TAB_UIR_RP <- ALL_TAB[,c("UIR_EE_S_RP_N",
                             "UIR_EE_T_RP_N",
                             "UIR_NE_S_RP_N",
                             "UIR_NE_T_RP_N",
                             "UIR_EN_S_RP_N",
                             "UIR_EN_T_RP_N",
                             "UIR_NN_S_RP_N",
                             "UIR_NN_T_RP_N"
  )]
  function.get_simple_twisted_boxplot(
    "Simple and twisted read pairs in undirected reference interactions",
    "Interaction number",
    undirected_ref_color,
    N_TAB_UIR_RP,
    YL)
  
  plot.new()
  
  # Differences between simple and twisted
  
    DIFF_VEC_N_DI <- function.get_simple_twisted_diff_boxplot(
    "Differences of directed simple and twisted interaction numbers",
    "Difference of interaction numbers",
    directed_color,
    N_TAB_DI,
    FALSE)
  
    DIFF_VEC_N_DI_RP <- function.get_simple_twisted_diff_boxplot(
    "Differences of simple and twisted read pairs in directed interactions",
    "Difference of interaction numbers",
    directed_color,
    N_TAB_DI_RP,
    FALSE)
  
    DIFF_VEC_N_UIR_RP <- function.get_simple_twisted_diff_boxplot(
    "Differences of simple and twisted read pairs in undirected reference interactions",
    "Difference of interaction numbers",
    undirected_ref_color,
    N_TAB_UIR_RP,
    FALSE)
  
    DIFF_VEC_N_UI_RP <- function.get_simple_twisted_diff_boxplot(
    "Differences of simple and twisted read pairs in undirected interactions",
    "Difference of interaction numbers",
    undirected_color,
    N_TAB_UI_RP,
    FALSE)
    
    # Plot differences between simple and twisted with comparable y-axes
    YMIN <- min(c(DIFF_VEC_N_DI_RP,DIFF_VEC_N_UIR_RP,DIFF_VEC_N_UI_RP))
    YMAX <- max(c(DIFF_VEC_N_DI_RP,DIFF_VEC_N_UIR_RP,DIFF_VEC_N_UI_RP))
    YL <- c(YMIN,YMAX)
    
    plot.new()
    
    DIFF_VEC_N_DI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs in directed interactions",
      "Difference of interaction numbers",
      directed_color,
      N_TAB_DI_RP,
      YL)
    
    DIFF_VEC_N_UIR_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs in undirected reference interactions",
      "Difference of interaction numbers",
      undirected_ref_color,
      N_TAB_UIR_RP,
      YL)
    
    DIFF_VEC_N_UI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs in undirected interactions",
      "Difference of interaction numbers",
      undirected_color,
      N_TAB_UI_RP,
      YL)
    
    mtext(paste("Interaction numbers",MM_TITLE_SUFFIX,sep=""), outer = TRUE, cex = 2)
    
  
dev.off()


# Median interaction and read pair distances
# ------------------------------------------

cairo_pdf(paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_median.pdf",sep=""), width=24, height=16)
par(mfrow=c(4,4), oma = c(0, 0, 3, 0))

# Directed interactions
MED_TAB_DI <- ALL_TAB[,c("DI_EE_S_MED",
                         "DI_EE_T_MED",
                         "DI_NE_S_MED",
                         "DI_NE_T_MED",
                         "DI_EN_S_MED",
                         "DI_EN_T_MED",
                         "DI_NN_S_MED",
                         "DI_NN_T_MED"
)]
MED_TAB_DI[,"DI_NE_S_MED"] = -MED_TAB_DI[,"DI_NE_S_MED"]
MED_TAB_DI[,"DI_NE_T_MED"] = -MED_TAB_DI[,"DI_NE_T_MED"]
function.get_simple_twisted_boxplot(
  "Directed simple and twisted interactions",
  "Median distance",
  directed_color,
  MED_TAB_DI,
  FALSE)

# Read pair numbers
# -----------------

# Directed interactions
MED_TAB_DI_RP <- ALL_TAB[,c("DI_EE_S_RP_MED",
                            "DI_EE_T_RP_MED",
                            "DI_NE_S_RP_MED",
                            "DI_NE_T_RP_MED",
                            "DI_EN_S_RP_MED",
                            "DI_EN_T_RP_MED",
                            "DI_NN_S_RP_MED",
                            "DI_NN_T_RP_MED"
)]
MED_TAB_DI_RP[,"DI_NE_S_RP_MED"] = -MED_TAB_DI_RP[,"DI_NE_S_RP_MED"]
MED_TAB_DI_RP[,"DI_NE_T_RP_MED"] = -MED_TAB_DI_RP[,"DI_NE_T_RP_MED"]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted directed interactions",
  "Median distance",
  directed_color,
  MED_TAB_DI_RP,
  FALSE)

# Undirected reference interactions
MED_TAB_UIR_RP <- ALL_TAB[,c("UIR_EE_S_RP_MED",
                             "UIR_EE_T_RP_MED",
                             "UIR_NE_S_RP_MED",
                             "UIR_NE_T_RP_MED",
                             "UIR_EN_S_RP_MED",
                             "UIR_EN_T_RP_MED",
                             "UIR_NN_S_RP_MED",
                             "UIR_NN_T_RP_MED"
)]
MED_TAB_UIR_RP[,"UIR_NE_S_RP_MED"] = -MED_TAB_UIR_RP[,"UIR_NE_S_RP_MED"]
MED_TAB_UIR_RP[,"UIR_NE_T_RP_MED"] = -MED_TAB_UIR_RP[,"UIR_NE_T_RP_MED"]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected reference interactions",
  "Median distance",
  undirected_ref_color,
  MED_TAB_UIR_RP,
  FALSE)

# Undirected interactions
MED_TAB_UI_RP <- ALL_TAB[,c("UI_EE_S_RP_MED",
                            "UI_EE_T_RP_MED",
                            "UI_NE_S_RP_MED",
                            "UI_NE_T_RP_MED",
                            "UI_EN_S_RP_MED",
                            "UI_EN_T_RP_MED",
                            "UI_NN_S_RP_MED",
                            "UI_NN_T_RP_MED"
)]
MED_TAB_UI_RP[,"UI_NE_S_RP_MED"] = -MED_TAB_UI_RP[,"UI_NE_S_RP_MED"]
MED_TAB_UI_RP[,"UI_NE_T_RP_MED"] = -MED_TAB_UI_RP[,"UI_NE_T_RP_MED"]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected interactions",
  "Median distance",
  undirected_color,
  MED_TAB_UI_RP,
  FALSE)

YL <- c(
  min(MED_TAB_DI_RP[,c("DI_NE_S_RP_MED","DI_NE_T_RP_MED","DI_EN_S_RP_MED","DI_EN_T_RP_MED")],
      MED_TAB_UIR_RP[,c("UIR_NE_S_RP_MED","UIR_NE_T_RP_MED","UIR_EN_S_RP_MED","UIR_EN_T_RP_MED")]),
  max(MED_TAB_DI_RP[,c("DI_NE_S_RP_MED","DI_NE_T_RP_MED","DI_EN_S_RP_MED","DI_EN_T_RP_MED")],
      MED_TAB_UIR_RP[,c("UIR_NE_S_RP_MED","UIR_NE_T_RP_MED","UIR_EN_S_RP_MED","UIR_EN_T_RP_MED")])
)

plot.new()

function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted directed interactions",
  "Median distance",
  directed_color,
  MED_TAB_DI_RP,
  YL)

function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected reference interactions",
  "Median distance",
  undirected_ref_color,
  MED_TAB_UIR_RP,
  YL)

plot.new()

# Differences between simple and twisted
DIFF_VEC_MED_DI <- function.get_simple_twisted_diff_boxplot(
  "Differences for directed simple and twisted interactions",
  "Difference of median distances",
  directed_color,
  MED_TAB_DI,
  FALSE)

DIFF_VEC_MED_DI_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in directed interactions",
  "Difference of median distances",
  directed_color,
  MED_TAB_DI_RP,
  FALSE)

DIFF_VEC_MED_UIR_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in undirected reference interactions",
  "Difference of median distances",
  undirected_ref_color,
  MED_TAB_UIR_RP,
  FALSE)

DIFF_VEC_MED_UI_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in undirected interactions",
  "Difference of median distances",
  undirected_color,
  MED_TAB_UI_RP,
  FALSE)

# Plot differences between simple and twisted with comparable y-axes
YMIN <- min(c(DIFF_VEC_MED_DI_RP,DIFF_VEC_MED_UIR_RP,DIFF_VEC_MED_UI_RP))
YMAX <- max(c(DIFF_VEC_MED_DI_RP,DIFF_VEC_MED_UIR_RP,DIFF_VEC_MED_UI_RP))
YL <- c(YMIN,YMAX)

plot.new()

DIFF_VEC_MED_DI_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in directed interactions",
  "Difference of median distances",
  directed_color,
  MED_TAB_DI_RP,
  YL)

DIFF_VEC_MED_UIR_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in undirected reference interactions",
  "Difference of median distances",
  undirected_ref_color,
  MED_TAB_UIR_RP,
  YL)

DIFF_VEC_MED_UI_RP <- function.get_simple_twisted_diff_boxplot(
  "Differences for simple and twisted read pairs in undirected interactions",
  "Difference of median distances",
  undirected_color,
  MED_TAB_UI_RP,
  YL)

  mtext(paste("Median distances",MM_TITLE_SUFFIX,sep=""), outer = TRUE, cex = 2)

dev.off()


# Interquartile ranges of interaction and read pair distances
# ------------------------------------------

cairo_pdf(paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_iqr.pdf",sep=""), width=24, height=16)
par(mfrow=c(4,4), oma = c(0, 0, 3, 0))

# Directed interactions
IQR_TAB_DI <- ALL_TAB[,c("DI_EE_S_IQR",
                         "DI_EE_T_IQR",
                         "DI_NE_S_IQR",
                         "DI_NE_T_IQR",
                         "DI_EN_S_IQR",
                         "DI_EN_T_IQR",
                         "DI_NN_S_IQR",
                         "DI_NN_T_IQR"
)]
function.get_simple_twisted_boxplot(
  "IQR of distances of directed simple and twisted interactions",
  "IQR of interaction distance",
  directed_color,
  IQR_TAB_DI,
  FALSE)

# Read pair numbers
# -----------------

# Directed interactions
IQR_TAB_DI_RP <- ALL_TAB[,c("DI_EE_S_RP_IQR",
                            "DI_EE_T_RP_IQR",
                            "DI_NE_S_RP_IQR",
                            "DI_NE_T_RP_IQR",
                            "DI_EN_S_RP_IQR",
                            "DI_EN_T_RP_IQR",
                            "DI_NN_S_RP_IQR",
                            "DI_NN_T_RP_IQR"
)]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted directed interactions",
  "IQR of interaction distance",
  directed_color,
  IQR_TAB_DI_RP,
  FALSE)

# Undirected reference interactions
IQR_TAB_UIR_RP <- ALL_TAB[,c("UIR_EE_S_RP_IQR",
                             "UIR_EE_T_RP_IQR",
                             "UIR_NE_S_RP_IQR",
                             "UIR_NE_T_RP_IQR",
                             "UIR_EN_S_RP_IQR",
                             "UIR_EN_T_RP_IQR",
                             "UIR_NN_S_RP_IQR",
                             "UIR_NN_T_RP_IQR"
)]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected reference interactions",
  "IQR of interaction distance",
  undirected_ref_color,
  IQR_TAB_UIR_RP,
  FALSE)

# Undirected interactions
IQR_TAB_UI_RP <- ALL_TAB[,c("UI_EE_S_RP_IQR",
                            "UI_EE_T_RP_IQR",
                            "UI_NE_S_RP_IQR",
                            "UI_NE_T_RP_IQR",
                            "UI_EN_S_RP_IQR",
                            "UI_EN_T_RP_IQR",
                            "UI_NN_S_RP_IQR",
                            "UI_NN_T_RP_IQR"
)]
function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected interactions",
  "IQR of interaction distance",
  undirected_color,
  IQR_TAB_UI_RP,
  FALSE)

YL <- c(
  min(IQR_TAB_DI_RP[,c("DI_NE_S_RP_IQR","DI_NE_T_RP_IQR","DI_EN_S_RP_IQR","DI_EN_T_RP_IQR")],
      IQR_TAB_UIR_RP[,c("UIR_NE_S_RP_IQR","UIR_NE_T_RP_IQR","UIR_EN_S_RP_IQR","UIR_EN_T_RP_IQR")]),
  max(IQR_TAB_DI_RP[,c("DI_NE_S_RP_IQR","DI_NE_T_RP_IQR","DI_EN_S_RP_IQR","DI_EN_T_RP_IQR")],
      IQR_TAB_UIR_RP[,c("UIR_NE_S_RP_IQR","UIR_NE_T_RP_IQR","UIR_EN_S_RP_IQR","UIR_EN_T_RP_IQR")])
)

plot.new()

function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted directed interactions",
  "IQR of interaction distance",
  directed_color,
  IQR_TAB_DI_RP,
  YL)

function.get_simple_twisted_boxplot(
  "Read pairs in simple and twisted undirected reference interactions",
  "IQR of interaction distance",
  undirected_ref_color,
  IQR_TAB_UIR_RP,
  YL)

plot.new()

# Differences between simple and twisted
DIFF_VEC_IQR_DI <- function.get_simple_twisted_diff_boxplot(
  "Directed simple and twisted interactions",
  "Difference of IQR of distances",
  directed_color,
  IQR_TAB_DI,
  FALSE)

DIFF_VEC_IQR_DI_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in directed interactions",
  "Difference of IQR of distances",
  directed_color,
  IQR_TAB_DI_RP,
  FALSE)

DIFF_VEC_IQR_UIR_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in undirected reference interactions",
  "Difference of IQR of distances",
  undirected_ref_color,
  IQR_TAB_UIR_RP,
  FALSE)

DIFF_VEC_IQR_UI_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in undirected interactions",
  "Difference of IQR of distances",
  undirected_color,
  IQR_TAB_UI_RP,
  FALSE)

# Plot differences between simple and twisted with comparable y-axes
YMIN <- min(c(DIFF_VEC_IQR_DI_RP,DIFF_VEC_IQR_UIR_RP,DIFF_VEC_IQR_UI_RP))
YMAX <- max(c(DIFF_VEC_IQR_DI_RP,DIFF_VEC_IQR_UIR_RP,DIFF_VEC_IQR_UI_RP))
YL <- c(YMIN,YMAX)

plot.new()

DIFF_VEC_IQR_DI_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in directed interactions",
  "Difference of IQR of distances",
  directed_color,
  IQR_TAB_DI_RP,
  YL)

DIFF_VEC_IQR_UIR_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in undirected reference interactions",
  "Difference of IQR of distances",
  undirected_ref_color,
  IQR_TAB_UIR_RP,
  YL)

DIFF_VEC_IQR_UI_RP <- function.get_simple_twisted_diff_boxplot(
  "Distances of simple and twisted read pairs in undirected interactions",
  "Difference of IQR of distances",
  undirected_color,
  IQR_TAB_UI_RP,
  YL)

  mtext(paste("Interquartile ranges",MM_TITLE_SUFFIX,sep=""), outer = TRUE, cex = 2)

dev.off()
