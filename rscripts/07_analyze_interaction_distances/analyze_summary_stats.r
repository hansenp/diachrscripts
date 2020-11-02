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


function.get_pdf_for_n_med_or_iqr <- function(
  PDF_NAME,
  MM_TITLE,
  YLAB_1,
  YLAB_2,
  DI_TAB,
  UIR_TAB,
  UI_TAB
)
{
  # Remove tag for measure ('_N','_MED' or '_IQR')
  colnames(DI_TAB) <- c("DI_EE","DI_NE","DI_EN","DI_NN")
  colnames(UIR_TAB) <- c("UIR_EE","UIR_NE","UIR_EN","UIR_NN")
  colnames(UI_TAB) <- c("UI_EE","UI_NE","UI_EN","UI_NN")
  
  # Init PDF file
  cairo_pdf(PDF_NAME, width=18, height=21)
  par(mfrow=c(7,4), oma = c(0, 0, 3, 0))
  
  # Create one boxplot for each interaction category
  boxplot(
    DI_TAB,
    main="Directed interactions",
    xlab="Enrichment status",
    ylab=YLAB_1,
    col=directed_color,
    names=c("EE", "NE", "EN", "NN")
  )
  boxplot(
    UIR_TAB,
    main="Directed interactions",
    xlab="Enrichment status",
    ylab=YLAB_1,
    col=undirected_ref_color,
    names=c("EE", "NE", "EN", "NN")
  )
  boxplot(
    UI_TAB,
    main="Directed interactions",
    xlab="Enrichment status",
    ylab=YLAB_1,
    col=undirected_color,
    names=c("EE", "NE", "EN", "NN")
  )
  
  plot.new()
  
  # Create one boxplot for all interaction categories
  ALL_TAB <- cbind(DI_TAB,UIR_TAB,UI_TAB)
  boxplot(
    ALL_TAB,
    main="Directed, undirected reference and undirected interactions",
    xlab="Enrichment status",
    ylab=YLAB_1,  
    col=c(directed_color,directed_color,directed_color,directed_color,
          undirected_ref_color,undirected_ref_color,undirected_ref_color,undirected_ref_color,
          undirected_color,undirected_color,undirected_color,undirected_color),
    names=c("EE", "NE", "EN", "NN","EE", "NE", "EN", "NN","EE", "NE", "EN", "NN")
  )
  
  # Create boxplot for interaction categories DI and UIR only
  TMP_TAB <- ALL_TAB[,c("DI_EE","UIR_EE","DI_NE","UIR_NE","DI_EN","UIR_EN","DI_NN","UIR_NN")]
  boxplot(
    TMP_TAB,
    main="Directed and undirected reference interactions",
    xlab="Enrichment status",
    ylab=YLAB_1,
    col=c(directed_color,undirected_ref_color),
    names=c("EE", "EE", "NE", "NE","EN", "EN", "NN", "NN")
  )
  
  plot.new()
  plot.new()
  
  # Compare DI with UIR using scatterplots
  
  YMIN <- min(ALL_TAB[,c("DI_EE","DI_NE","DI_EN","DI_NN","UIR_EE","UIR_NE","UIR_EN","UIR_NN")])
  YMAX <- max(ALL_TAB[,c("DI_EE","DI_NE","DI_EN","DI_NN","UIR_EE","UIR_NE","UIR_EN","UIR_NN")])
  XMIN<-YMIN
  XMAX<-YMAX
  
  plot(
    ALL_TAB[,"DI_EE"]~ALL_TAB[,"UIR_EE"],
    main="Directed interactions - DI vs. UIR - EE",
    col=mixed_color_di_uir,
    xlim=c(XMIN,XMAX),
    ylim=c(YMIN,YMAX),
    xlab="UIR - EE",
    ylab="DI - EE"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"DI_EE"]~ALL_TAB[,"UIR_EE"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
  
  plot(
    ALL_TAB[,"DI_NE"]~ALL_TAB[,"UIR_NE"],
    main="Directed interactions - DI vs. UIR - NE",
    col=mixed_color_di_uir,
    xlim=c(XMIN,XMAX),
    ylim=c(YMIN,YMAX),
    xlab="UIR - NE",
    ylab="DI - NE"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"DI_NE"]~ALL_TAB[,"UIR_NE"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
  
  plot(
    ALL_TAB[,"DI_EN"]~ALL_TAB[,"UIR_EN"],
    main="Directed interactions - DI vs. UIR - EN",
    col=mixed_color_di_uir,
    xlim=c(XMIN,XMAX),
    ylim=c(YMIN,YMAX),
    xlab="UIR - EN",
    ylab="DI - EN"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"DI_EN"]~ALL_TAB[,"UIR_EN"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
  
  plot(
    ALL_TAB[,"DI_NN"]~ALL_TAB[,"UIR_NN"],
    main="Directed interactions - DI vs. UIR - NN",
    col=mixed_color_di_uir,
    xlim=c(XMIN,XMAX),
    ylim=c(YMIN,YMAX),
    xlab="UIR - NN",
    ylab="DI - NN"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"DI_NN"]~ALL_TAB[,"UIR_NN"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
  
  
  # Get differences of DI and UIR in different enrichment categories
  DIFF_DI_UIR_EE <- ALL_TAB[,"DI_EE"]-ALL_TAB[,"UIR_EE"]
  DIFF_DI_UIR_NE <- ALL_TAB[,"DI_NE"]-ALL_TAB[,"UIR_NE"]
  DIFF_DI_UIR_EN <- ALL_TAB[,"DI_EN"]-ALL_TAB[,"UIR_EN"]
  DIFF_DI_UIR_NN <- ALL_TAB[,"DI_NN"]-ALL_TAB[,"UIR_NN"]
  
  # Get minimum and maximum values required for postioning of P-values
  y_min <- min(DIFF_DI_UIR_EE,DIFF_DI_UIR_NE,DIFF_DI_UIR_EN,DIFF_DI_UIR_NN)
  y_max <- max(DIFF_DI_UIR_EE,DIFF_DI_UIR_NE,DIFF_DI_UIR_EN,DIFF_DI_UIR_NN)
  y_range <- abs(y_min-y_max)
  YMIN <- y_min
  YMAX <- y_max + y_range/10
  
  # Create boxplot for differences
  boxplot(
    cbind(DIFF_DI_UIR_EE,DIFF_DI_UIR_NE,DIFF_DI_UIR_EN,DIFF_DI_UIR_NN),
    main="Differences of DI and UIR interactions",
    xlab="Interaction category",
    ylab=YLAB_2, 
    col=c(mixed_color_di_uir,mixed_color_di_uir,mixed_color_di_uir,mixed_color_di_uir),
    names=c("EE", "NE", "EN", "NN"),
    ylim=c(YMIN,YMAX)
  )
  abline(h=0,col="grey")
  
  # Perform t-tests on differences
  DIFF_DI_UIR_EE_TT <- t.test(DIFF_DI_UIR_EE)
  DIFF_DI_UIR_NE_TT <- t.test(DIFF_DI_UIR_NE)
  DIFF_DI_UIR_EN_TT <- t.test(DIFF_DI_UIR_EN)
  DIFF_DI_UIR_NN_TT <- t.test(DIFF_DI_UIR_NN)
  
  # Add P-values to boxplot
  v_space <- y_range/20
  text(c(1:4),
       c(max(DIFF_DI_UIR_EE)+v_space,
         max(DIFF_DI_UIR_NE)+v_space,
         max(DIFF_DI_UIR_EN)+v_space,
         max(DIFF_DI_UIR_NN)+v_space),
       c(paste("p = ", formatC(DIFF_DI_UIR_EE_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(DIFF_DI_UIR_NE_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(DIFF_DI_UIR_EN_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(DIFF_DI_UIR_NN_TT$p.value, format = "e", digits = 2), sep="")
       )
  )
  
  plot.new()
  plot.new()
  plot.new()
  
  # Create boxplot for NE and EN of directed and undirected reference interactions
  boxplot(
    ALL_TAB[,c("DI_NE","DI_EN")],
    main="Directed interactions - NE and EN",
    xlab="Enrichment status",
    ylab=YLAB_1,  
    col=c(directed_color,directed_color),
    names=c("NE", "EN")
  )
  boxplot(
    ALL_TAB[,c("UIR_NE","UIR_EN")],
    main="Undirected reference interactions - NE and EN",
    xlab="Enrichment status",
    ylab=YLAB_1,  
    col=c(undirected_ref_color,undirected_ref_color),
    names=c("NE", "EN")
  )
  boxplot(
    ALL_TAB[,c("UI_NE","UI_EN")],
    main="Undirected interactions - NE and EN",
    xlab="Enrichment status",
    ylab=YLAB_1,  
    col=c(undirected_color,undirected_color),
    names=c("NE", "EN")
  )
  
  plot.new()
  
  plot(
    ALL_TAB[,"DI_NE"]~ALL_TAB[,"DI_EN"],
    main="Directed interactions - NE vs. EN",
    col=directed_color,
    xlim=c(min(ALL_TAB[,"DI_NE"],ALL_TAB[,"DI_EN"]),max(ALL_TAB[,"DI_NE"],ALL_TAB[,"DI_EN"])),
    ylim=c(min(ALL_TAB[,"DI_NE"],ALL_TAB[,"DI_EN"]),max(ALL_TAB[,"DI_NE"],ALL_TAB[,"DI_EN"])),
    xlab="EN",
    ylab="NE"
    )
  abline(0,1,col="gray")
  text(ALL_TAB[,"DI_NE"]~ALL_TAB[,"DI_EN"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)
  
  plot(
    ALL_TAB[,"UIR_NE"]~ALL_TAB[,"UIR_EN"],
    main="Undirected reference interactions - NE vs. EN",
    col=undirected_ref_color,
    xlim=c(min(ALL_TAB[,"UIR_NE"],ALL_TAB[,"UIR_EN"]),max(ALL_TAB[,"UIR_NE"],ALL_TAB[,"UIR_EN"])),
    ylim=c(min(ALL_TAB[,"UIR_NE"],ALL_TAB[,"UIR_EN"]),max(ALL_TAB[,"UIR_NE"],ALL_TAB[,"UIR_EN"])),
    xlab="EN",
    ylab="NE"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"UIR_NE"]~ALL_TAB[,"UIR_EN"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)

  plot(
    ALL_TAB[,"UI_NE"]~ALL_TAB[,"UI_EN"],
    main="Undirected interactions - NE vs. EN",
    col=undirected_color,
    xlim=c(min(ALL_TAB[,"UI_NE"],ALL_TAB[,"UI_EN"]),max(ALL_TAB[,"UI_NE"],ALL_TAB[,"UI_EN"])),
    ylim=c(min(ALL_TAB[,"UI_NE"],ALL_TAB[,"UI_EN"]),max(ALL_TAB[,"UI_NE"],ALL_TAB[,"UI_EN"])),
    xlab="EN",
    ylab="NE"
  )
  abline(0,1,col="gray")
  text(ALL_TAB[,"UI_NE"]~ALL_TAB[,"UI_EN"], labels=rownames(ALL_TAB),data=ALL_TAB, cex=0.3)

  plot.new()

  # Get differences of NE and EN for different interaction categories
  DIFF_NE_EN_DI <- ALL_TAB[,"DI_NE"]-ALL_TAB[,"DI_EN"]
  DIFF_NE_EN_UIR <- ALL_TAB[,"UIR_NE"]-ALL_TAB[,"UIR_EN"]
  DIFF_NE_EN_UI <- ALL_TAB[,"UI_NE"]-ALL_TAB[,"UI_EN"]
  
  # Get minimum and maximum values required for postioning of P-values
  y_min <- min(DIFF_NE_EN_DI,DIFF_NE_EN_UIR,DIFF_NE_EN_UI)
  y_max <- max(DIFF_NE_EN_DI,DIFF_NE_EN_UIR,DIFF_NE_EN_UI)
  y_range <- abs(y_min-y_max)
  YMIN <- y_min
  YMAX <- y_max + y_range/10
  
  # Create boxplot for differences
  boxplot(
    cbind(DIFF_NE_EN_DI,DIFF_NE_EN_UIR,DIFF_NE_EN_UI),
    main="Differences of NE and EN interactions",
    xlab="Interaction category",
    ylab=YLAB_2, 
    col=c(directed_color,undirected_ref_color,undirected_color),
    names=c("DI", "UIR", "UI"),
    ylim=c(YMIN,YMAX)
  )
  abline(h=0,col="grey")
  
  # Perform t-tests on differences
  DIFF_NE_EN_DI_TT <- t.test(DIFF_NE_EN_DI)
  DIFF_NE_EN_UIR_TT <- t.test(DIFF_NE_EN_UIR)
  DIFF_NE_EN_UI_TT <- t.test(DIFF_NE_EN_UI)
  
  # Add P-values to boxplot
  v_space <- y_range/20
  text(c(1:3),
       c(max(DIFF_NE_EN_DI)+v_space,
         max(DIFF_NE_EN_UIR)+v_space,
         max(DIFF_NE_EN_UI)+v_space),
       c(paste("p = ", formatC(DIFF_NE_EN_DI_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(DIFF_NE_EN_UIR_TT$p.value, format = "e", digits = 2), sep=""),
         paste("p = ", formatC(DIFF_NE_EN_UI_TT$p.value, format = "e", digits = 2), sep="")
       )
  )
  
  # Add title above all plots
  mtext(MM_TITLE, outer = TRUE, cex = 2)
  
  # Close PDF file
  dev.off()
}

N_TAB_DI <- ALL_TAB[,c("DI_EE_N","DI_NE_N","DI_EN_N","DI_NN_N")]
N_TAB_UIR <- ALL_TAB[,c("UIR_EE_N","UIR_NE_N","UIR_EN_N","UIR_NN_N")]
N_TAB_UI <- ALL_TAB[,c("UI_EE_N","UI_NE_N","UI_EN_N","UI_NN_N")]
function.get_pdf_for_n_med_or_iqr(
  paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_n.pdf",sep=""),
  paste("Interaction numbers",MM_TITLE_SUFFIX,sep=""),
  "Interaction number",
  "Difference of interaction numbers",
  N_TAB_DI,
  N_TAB_UIR,
  N_TAB_UI  
)

MED_TAB_DI <- ALL_TAB[,c("DI_EE_MED","DI_NE_MED","DI_EN_MED","DI_NN_MED")]
MED_TAB_UIR <- ALL_TAB[,c("UIR_EE_MED","UIR_NE_MED","UIR_EN_MED","UIR_NN_MED")]
MED_TAB_UI <- ALL_TAB[,c("UI_EE_MED","UI_NE_MED","UI_EN_MED","UI_NN_MED")]
function.get_pdf_for_n_med_or_iqr(
  paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_median.pdf",sep=""),
  paste("Median distances",MM_TITLE_SUFFIX,sep=""),
  "Median distance",
  "Difference of median distances",
  MED_TAB_DI,
  MED_TAB_UIR,
  MED_TAB_UI  
)

IQR_TAB_DI <- ALL_TAB[,c("DI_EE_IQR","DI_NE_IQR","DI_EN_IQR","DI_NN_IQR")]
IQR_TAB_UIR <- ALL_TAB[,c("UIR_EE_IQR","UIR_NE_IQR","UIR_EN_IQR","UIR_NN_IQR")]
IQR_TAB_UI <- ALL_TAB[,c("UI_EE_IQR","UI_NE_IQR","UI_EN_IQR","UI_NN_IQR")]
function.get_pdf_for_n_med_or_iqr(
  paste(OUT_DIR,OUT_PREFIX,"interaction_distance_summary_stats_iqr.pdf",sep=""),
  paste("Interquartile ranges",MM_TITLE_SUFFIX,sep=""),
  "Interquartile range",
  "Difference of interquartile ranges",
  IQR_TAB_DI,
  IQR_TAB_UIR,
  IQR_TAB_UI  
)
