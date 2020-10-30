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

rownames(ALL_TAB) <- c("MK",
                       "ERY",
                       "NEU",
                       "MON",
                       "MAC_M0",
                       "MAC_M1",
                       "MAC_M2",
                       "EP",
                       "NB",
                       "TB",
                       "FOET",
                       "NCD4",
                       "TCD4",
                       "NACD4",
                       "ACD4",
                       "NCD8",
                       "TCD8")

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
  return(c(ST_DIFF_EE,ST_DIFF_NE,ST_DIFF_EN,ST_DIFF_NN))

}

function.get_pdf_for_n_med_or_iqr_st <- function(
  PDF_NAME,
  MM_TITLE,
  YLAB_1,
  YLAB_2,
  YLAB_3,
  YLAB_4,
  DI_TAB,
  DI_RP_TAB,
  UIR_RP_TAB,
  UI_RP_TAB
)
{
  cairo_pdf(PDF_NAME, width=18, height=30)
  par(mfrow=c(10,4), oma = c(0, 0, 3, 0))
  
    # ------------------------------
  
    function.get_simple_twisted_boxplot(
      "Simple and twisted DI",
      YLAB_1,
      directed_color,
      DI_TAB,
      FALSE)
  
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from DI",
      YLAB_2,
      directed_color,
      DI_RP_TAB,
      FALSE)
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from UIR",
      YLAB_2,
      undirected_ref_color,
      UIR_RP_TAB,
      FALSE)
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from UI",
      YLAB_2,
      undirected_color,
      UI_RP_TAB,
      FALSE)
    
    # ------------------------------
    
    y_min <- min(DI_TAB,DI_RP_TAB,UIR_RP_TAB,UI_RP_TAB)
    y_max <- max(DI_TAB,DI_RP_TAB,UIR_RP_TAB,UI_RP_TAB)
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted DI",
      YLAB_1,
      directed_color,
      DI_TAB,
      c(y_min,y_max))
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from DI",
      YLAB_2,
      directed_color,
      DI_RP_TAB,
      c(y_min,y_max))
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from UIR",
      YLAB_2,
      undirected_ref_color,
      UIR_RP_TAB,
      c(y_min,y_max))
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from UI",
      YLAB_2,
      undirected_color,
      UI_RP_TAB,
      c(y_min,y_max))
    
    # ------------------------------
    
    plot.new()
    
    y_min <- min(DI_RP_TAB,UIR_RP_TAB)
    y_max <- max(DI_RP_TAB,UIR_RP_TAB)
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from DI",
      YLAB_2,
      directed_color,
      DI_RP_TAB,
      c(y_min,y_max))
    
    function.get_simple_twisted_boxplot(
      "Simple and twisted read pairs from UIR",
      YLAB_2,
      undirected_ref_color,
      UIR_RP_TAB,
      c(y_min,y_max))

    plot.new()
    
    # ------------------------------
    
    YMIN <- min(DI_TAB)
    YMAX <- max(DI_TAB)
    XMIN<-YMIN
    XMAX<-YMAX
    
    plot(
      DI_TAB[,1]~DI_TAB[,2],
      main="Simple vs. twisted DI - EE",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_TAB[,1]~DI_TAB[,2], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_TAB[,3]~DI_TAB[,4],
      main="Simple vs. twisted DI - NE",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_TAB[,3]~DI_TAB[,4], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_TAB[,5]~DI_TAB[,6],
      main="Simple vs. twisted DI - EN",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_TAB[,5]~DI_TAB[,6], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_TAB[,7]~DI_TAB[,8],
      main="Simple vs. twisted DI - NN",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_TAB[,7]~DI_TAB[,8], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    # ------------------------------
    
    YMIN <- min(DI_RP_TAB)
    YMAX <- max(DI_RP_TAB)
    XMIN<-YMIN
    XMAX<-YMAX
    
    plot(
      DI_RP_TAB[,1]~DI_RP_TAB[,2],
      main="Simple vs. twisted read pairs from DI - EE",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_RP_TAB[,1]~DI_RP_TAB[,2], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_RP_TAB[,3]~DI_RP_TAB[,4],
      main="Simple vs. twisted read pairs from DI - NE",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_RP_TAB[,3]~DI_RP_TAB[,4], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_RP_TAB[,5]~DI_RP_TAB[,6],
      main="Simple vs. twisted read pairs from DI - EN",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_RP_TAB[,5]~DI_RP_TAB[,6], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      DI_RP_TAB[,7]~DI_RP_TAB[,8],
      main="Simple vs. twisted read pairs from DI - NN",
      col=directed_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(DI_RP_TAB[,7]~DI_RP_TAB[,8], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    # ------------------------------
    
    YMIN <- min(UIR_RP_TAB)
    YMAX <- max(UIR_RP_TAB)
    XMIN<-YMIN
    XMAX<-YMAX
    
    plot(
      UIR_RP_TAB[,1]~UIR_RP_TAB[,2],
      main="Simple vs. twisted read pairs from UIR - EE",
      col=undirected_ref_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UIR_RP_TAB[,1]~UIR_RP_TAB[,2], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UIR_RP_TAB[,3]~UIR_RP_TAB[,4],
      main="Simple vs. twisted read pairs from UIR - NE",
      col=undirected_ref_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UIR_RP_TAB[,3]~UIR_RP_TAB[,4], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UIR_RP_TAB[,5]~UIR_RP_TAB[,6],
      main="Simple vs. twisted read pairs from UIR - EN",
      col=undirected_ref_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UIR_RP_TAB[,5]~UIR_RP_TAB[,6], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UIR_RP_TAB[,7]~UIR_RP_TAB[,8],
      main="Simple vs. twisted read pairs from UIR - NN",
      col=undirected_ref_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UIR_RP_TAB[,7]~UIR_RP_TAB[,8], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    # ------------------------------
    
    YMIN <- min(UI_RP_TAB)
    YMAX <- max(UI_RP_TAB)
    XMIN<-YMIN
    XMAX<-YMAX
    
    plot(
      UI_RP_TAB[,1]~UI_RP_TAB[,2],
      main="Simple vs. twisted read pairs from UIR - EE",
      col=undirected_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UI_RP_TAB[,1]~UI_RP_TAB[,2], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UI_RP_TAB[,3]~UI_RP_TAB[,4],
      main="Simple vs. twisted read pairs from UIR - NE",
      col=undirected_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UI_RP_TAB[,3]~UI_RP_TAB[,4], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UI_RP_TAB[,5]~UI_RP_TAB[,6],
      main="Simple vs. twisted read pairs from UIR - EN",
      col=undirected_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UI_RP_TAB[,5]~UI_RP_TAB[,6], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    plot(
      UI_RP_TAB[,7]~UI_RP_TAB[,8],
      main="Simple vs. twisted read pairs from UIR - NN",
      col=undirected_color,
      xlim=c(XMIN,XMAX),
      ylim=c(YMIN,YMAX),
      xlab="Twisted",
      ylab="Simple"
    )
    abline(0,1,col="gray")
    text(UI_RP_TAB[,7]~UI_RP_TAB[,8], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
    
    # ------------------------------
    
    DIFF_VEC_DI <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted DI",
      YLAB_3,
      directed_color,
      DI_TAB,
      FALSE)
    
    DIFF_VEC_DI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      directed_color,
      DI_RP_TAB,
      FALSE)
    
    DIFF_VEC_UIR_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_ref_color,
        UIR_RP_TAB,
      FALSE)
    
    DIFF_VEC_UI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_color,
      UI_RP_TAB,
      FALSE)

    # ------------------------------
    
    YMIN <- min(DIFF_VEC_DI,DIFF_VEC_DI_RP,DIFF_VEC_UIR_RP,DIFF_VEC_UI_RP)
    YMAX <- max(DIFF_VEC_DI,DIFF_VEC_DI_RP,DIFF_VEC_UIR_RP,DIFF_VEC_UI_RP)
    XMIN<-YMIN
    XMAX<-YMAX
    
    DIFF_VEC_DI <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted DI",
      YLAB_3,
      directed_color,
      DI_TAB,
      c(YMIN,YMAX))
    
    DIFF_VEC_DI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      directed_color,
      DI_RP_TAB,
      c(YMIN,YMAX))
    
    DIFF_VEC_UIR_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_ref_color,
      UIR_RP_TAB,
      c(YMIN,YMAX))
    
    DIFF_VEC_UI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_color,
      UI_RP_TAB,
      c(YMIN,YMAX))
    
    # ------------------------------
    
    YMIN <- min(DIFF_VEC_DI_RP,DIFF_VEC_UIR_RP,DIFF_VEC_UI_RP)
    YMAX <- max(DIFF_VEC_DI_RP,DIFF_VEC_UIR_RP,DIFF_VEC_UI_RP)
    XMIN<-YMIN
    XMAX<-YMAX
    
    plot.new()
    
    DIFF_VEC_DI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      directed_color,
      DI_RP_TAB,
      c(YMIN,YMAX))
    
    DIFF_VEC_UIR_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_ref_color,
      UIR_RP_TAB,
      c(YMIN,YMAX))
    
    DIFF_VEC_UI_RP <- function.get_simple_twisted_diff_boxplot(
      "Differences of simple and twisted read pairs from DI",
      YLAB_4,
      undirected_color,
      UI_RP_TAB,
      c(YMIN,YMAX))
    
    
    mtext(MM_TITLE, outer=TRUE, cex=2)
  
  dev.off()
}

N_TAB_DI <- ALL_TAB[,c(
  "DI_EE_S_N","DI_EE_T_N",
  "DI_NE_S_N","DI_NE_T_N",
  "DI_EN_S_N","DI_EN_T_N",
  "DI_NN_S_N","DI_NN_T_N")]
N_TAB_DI_RP <- ALL_TAB[,c(
  "DI_EE_S_RP_N","DI_EE_T_RP_N",
  "DI_NE_S_RP_N","DI_NE_T_RP_N",
  "DI_EN_S_RP_N","DI_EN_T_RP_N",
  "DI_NN_S_RP_N","DI_NN_T_RP_N")]
N_TAB_UIR_RP <- ALL_TAB[,c(
  "UIR_EE_S_RP_N","UIR_EE_T_RP_N",
  "UIR_NE_S_RP_N","UIR_NE_T_RP_N",
  "UIR_EN_S_RP_N","UIR_EN_T_RP_N",
  "UIR_NN_S_RP_N","UIR_NN_T_RP_N")]
N_TAB_UI_RP <- ALL_TAB[,c(
  "UI_EE_S_RP_N","UI_EE_T_RP_N",
  "UI_NE_S_RP_N","UI_NE_T_RP_N",
  "UI_EN_S_RP_N","UI_EN_T_RP_N",
  "UI_NN_S_RP_N","UI_NN_T_RP_N")]

function.get_pdf_for_n_med_or_iqr_st(
  PDF_NAME=paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_n.pdf",sep=""),
  MM_TITLE=paste("Interaction and read pair numbers",MM_TITLE_SUFFIX,sep=""),
  YLAB_1="Interaction number",
  YLAB_2="Read pair number",
  YLAB_3="Difference of interaction numbers",
  YLAB_4="Difference of read pair numbers",
  DI_TAB=N_TAB_DI,
  DI_RP_TAB=N_TAB_DI_RP,
  UIR_RP_TAB=N_TAB_UIR_RP,
  UI_RP_TAB=N_TAB_UI_RP
)

MED_TAB_DI <- ALL_TAB[,c(
  "DI_EE_S_MED","DI_EE_T_MED",
  "DI_NE_S_MED","DI_NE_T_MED",
  "DI_EN_S_MED","DI_EN_T_MED",
  "DI_NN_S_MED","DI_NN_T_MED")]
MED_TAB_DI[,"DI_NE_S_MED"]<--MED_TAB_DI[,"DI_NE_S_MED"]
MED_TAB_DI[,"DI_NE_T_MED"]<--MED_TAB_DI[,"DI_NE_T_MED"]

MED_TAB_DI_RP <- ALL_TAB[,c(
  "DI_EE_S_RP_MED","DI_EE_T_RP_MED",
  "DI_NE_S_RP_MED","DI_NE_T_RP_MED",
  "DI_EN_S_RP_MED","DI_EN_T_RP_MED",
  "DI_NN_S_RP_MED","DI_NN_T_RP_MED")]
MED_TAB_DI_RP[,"DI_NE_S_RP_MED"]<--MED_TAB_DI_RP[,"DI_NE_S_RP_MED"]
MED_TAB_DI_RP[,"DI_NE_T_RP_MED"]<--MED_TAB_DI_RP[,"DI_NE_T_RP_MED"]

MED_TAB_UIR_RP <- ALL_TAB[,c(
  "UIR_EE_S_RP_MED","UIR_EE_T_RP_MED",
  "UIR_NE_S_RP_MED","UIR_NE_T_RP_MED",
  "UIR_EN_S_RP_MED","UIR_EN_T_RP_MED",
  "UIR_NN_S_RP_MED","UIR_NN_T_RP_MED")]
MED_TAB_UIR_RP[,"UIR_NE_S_RP_MED"]<--MED_TAB_UIR_RP[,"UIR_NE_S_RP_MED"]
MED_TAB_UIR_RP[,"UIR_NE_T_RP_MED"]<--MED_TAB_UIR_RP[,"UIR_NE_T_RP_MED"]

MED_TAB_UI_RP <- ALL_TAB[,c(
  "UI_EE_S_RP_MED","UI_EE_T_RP_MED",
  "UI_NE_S_RP_MED","UI_NE_T_RP_MED",
  "UI_EN_S_RP_MED","UI_EN_T_RP_MED",
  "UI_NN_S_RP_MED","UI_NN_T_RP_MED")]
MED_TAB_UI_RP[,"UI_NE_S_RP_MED"]<--MED_TAB_UI_RP[,"UI_NE_S_RP_MED"]
MED_TAB_UI_RP[,"UI_NE_T_RP_MED"]<--MED_TAB_UI_RP[,"UI_NE_T_RP_MED"]

function.get_pdf_for_n_med_or_iqr_st(
  PDF_NAME=paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_median.pdf",sep=""),
  MM_TITLE=paste("Median interaction and read pair distances",MM_TITLE_SUFFIX,sep=""),
  YLAB_1="Interaction distance",
  YLAB_2="Read pair distance",
  YLAB_3="Difference of median interaction distances",
  YLAB_4="Difference of median read pair distances",
  DI_TAB=MED_TAB_DI,
  DI_RP_TAB=MED_TAB_DI_RP,
  UIR_RP_TAB=MED_TAB_UIR_RP,
  UI_RP_TAB=MED_TAB_UI_RP
)


IQR_TAB_DI <- ALL_TAB[,c(
  "DI_EE_S_IQR","DI_EE_T_IQR",
  "DI_NE_S_IQR","DI_NE_T_IQR",
  "DI_EN_S_IQR","DI_EN_T_IQR",
  "DI_NN_S_IQR","DI_NN_T_IQR")]
IQR_TAB_DI_RP <- ALL_TAB[,c(
  "DI_EE_S_RP_IQR","DI_EE_T_RP_IQR",
  "DI_NE_S_RP_IQR","DI_NE_T_RP_IQR",
  "DI_EN_S_RP_IQR","DI_EN_T_RP_IQR",
  "DI_NN_S_RP_IQR","DI_NN_T_RP_IQR")]
IQR_TAB_UIR_RP <- ALL_TAB[,c(
  "UIR_EE_S_RP_IQR","UIR_EE_T_RP_IQR",
  "UIR_NE_S_RP_IQR","UIR_NE_T_RP_IQR",
  "UIR_EN_S_RP_IQR","UIR_EN_T_RP_IQR",
  "UIR_NN_S_RP_IQR","UIR_NN_T_RP_IQR")]
IQR_TAB_UI_RP <- ALL_TAB[,c(
  "UI_EE_S_RP_IQR","UI_EE_T_RP_IQR",
  "UI_NE_S_RP_IQR","UI_NE_T_RP_IQR",
  "UI_EN_S_RP_IQR","UI_EN_T_RP_IQR",
  "UI_NN_S_RP_IQR","UI_NN_T_RP_IQR")]

function.get_pdf_for_n_med_or_iqr_st(
  PDF_NAME=paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_st_iqr.pdf",sep=""),
  MM_TITLE=paste("Interquartile ranges (IQR) of interaction and read pair distances",MM_TITLE_SUFFIX,sep=""),
  YLAB_1="Interaction distance",
  YLAB_2="Read pair distance",
  YLAB_3="Difference of IQR of interaction distances",
  YLAB_4="Difference of IQR of read pair distances",
  DI_TAB=IQR_TAB_DI,
  DI_RP_TAB=IQR_TAB_DI_RP,
  UIR_RP_TAB=IQR_TAB_UIR_RP,
  UI_RP_TAB=IQR_TAB_UI_RP
)
