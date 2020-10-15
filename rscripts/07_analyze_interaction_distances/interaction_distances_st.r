#!/usr/bin/env Rscript
source("rscripts/07_analyze_interaction_distances/interaction_distances_lib.r")

# Analyze interaction distances with respect to simple and twisted orientation.
# At the level of interaction this can be done for directed interactions only.
# At the level of read pairs this can be done also for undirected interactions.

# Read commannd line arguments
args = commandArgs(trailingOnly=TRUE)

for(arg in args){
  print(arg)
}

OUT_DIR <- args[1]
OUT_PREFIX <- args[2]

DI_EE_S_DIST <- t(read.table(args[3]))
DI_NE_S_DIST <- t(-read.table(args[4]))
DI_EN_S_DIST <- t(read.table(args[5]))
DI_NN_S_DIST <- t(read.table(args[6]))

DI_EE_T_DIST <- t(read.table(args[7]))
DI_NE_T_DIST <- t(-read.table(args[8]))
DI_EN_T_DIST <- t(read.table(args[9]))
DI_NN_T_DIST <- t(read.table(args[10]))

DI_EE_S_RP_DIST <- t(read.table(args[11]))
DI_NE_S_RP_DIST <- t(-read.table(args[12]))
DI_EN_S_RP_DIST <- t(read.table(args[13]))
DI_NN_S_RP_DIST <- t(read.table(args[14]))

DI_EE_T_RP_DIST <- t(read.table(args[15]))
DI_NE_T_RP_DIST <- t(-read.table(args[16]))
DI_EN_T_RP_DIST <- t(read.table(args[17]))
DI_NN_T_RP_DIST <- t(read.table(args[18]))

UIR_EE_S_RP_DIST <- t(read.table(args[19]))
UIR_NE_S_RP_DIST <- t(-read.table(args[20]))
UIR_EN_S_RP_DIST <- t(read.table(args[21]))
UIR_NN_S_RP_DIST <- t(read.table(args[22]))

UIR_EE_T_RP_DIST <- t(read.table(args[23]))
UIR_NE_T_RP_DIST <- t(-read.table(args[24]))
UIR_EN_T_RP_DIST <- t(read.table(args[25]))
UIR_NN_T_RP_DIST <- t(read.table(args[26]))

UI_EE_S_RP_DIST <- t(read.table(args[27]))
UI_NE_S_RP_DIST <- t(-read.table(args[28]))
UI_EN_S_RP_DIST <- t(read.table(args[29]))
UI_NN_S_RP_DIST <- t(read.table(args[30]))

UI_EE_T_RP_DIST <- t(read.table(args[31]))
UI_NE_T_RP_DIST <- t(-read.table(args[32]))
UI_EN_T_RP_DIST <- t(read.table(args[33]))
UI_NN_T_RP_DIST <- t(read.table(args[34]))

function.plot_histograms_i_ee_ne_en_nn <- function(
  PDF_FILE_NAME,
  PDF_WIDTH,
  PDF_HEIGHT,
  MAIN_MAIN_TITLE,
  BREAKS,
  XMAX,
  BAR_COLOR,
  EE_S,
  NE_S,
  EN_S,
  NN_S,
  EE_T,
  NE_T,
  EN_T,
  NN_T 
) 
{
  # Get legend vectors with interaction numbers, median distances and interquartile ranges
  legend_vec_ee_s <- function.get_legend_vector(EE_S)
  legend_vec_ne_s <- function.get_legend_vector(NE_S)
  legend_vec_en_s <- function.get_legend_vector(EN_S)
  legend_vec_nn_s <- function.get_legend_vector(NN_S)
  
  legend_vec_ee_t <- function.get_legend_vector(EE_T)
  legend_vec_ne_t <- function.get_legend_vector(NE_T)
  legend_vec_en_t <- function.get_legend_vector(EN_T)
  legend_vec_nn_t <- function.get_legend_vector(NN_T)
  
  # Get common maximum for y-axes
  YMAX <- function.get_common_ymax_for_histograms(
    list(
      EE_S,
      NE_S,
      EN_S,
      NN_S,
      EE_T,
      NE_T,
      EN_T,
      NN_T     
    ),
    BREAKS)
  
  # Plot six histograms to PDF file
  cairo_pdf(PDF_FILE_NAME, width=PDF_WIDTH, height=PDF_HEIGHT)
  
  par(mfrow=c(4,5), oma = c(0, 0, 2, 0))
  
  # Simple
  hist(EE_S,
       main="Simple - EE",
       col=BAR_COLOR,
       border=simple_color, 
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_ee_s)
  
  hist(NE_S,
       main="Simple - NE (left)",
       col=BAR_COLOR,
       border=simple_color, 
       xlab="Interaction distance",
       xlim=c(-XMAX,0), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topleft", cex=0.8, bty="n",
         legend=legend_vec_ne_s)
  
  hist(EN_S,
       main="Simple - EN (right)",
       col=BAR_COLOR,
       border=simple_color, 
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_en_s)
  
  hist(NN_S,
       main="Simple - NN",
       col=BAR_COLOR,
       border=simple_color,
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_nn_s)
  
  # Get density differences of NE and EN for simple and twisted
  density_diff_list_ne_en_s <- function.get_density_diff(-NE_S,EN_S, BREAKS)
  density_diff_list_ne_en_t <- function.get_density_diff(-NE_T,EN_T, BREAKS)
  
  # Get common ylims of density differences for DI, UIR and UI
  MINY <- min(density_diff_list_ne_en_s$density_diff,
              density_diff_list_ne_en_t$density_diff
  )
  MAXY <- max(density_diff_list_ne_en_s$density_diff,
              density_diff_list_ne_en_t$density_diff
  )

  # Plot differences of density differences of NE and EN within simple
  plot(density_diff_list_ne_en_s$bin_centers,
       density_diff_list_ne_en_s$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - (NE-EN)",
       xlab="Interaction distance",
       ylab="NE minus EN",
       pch=21,
       col=BAR_COLOR,
       bg=simple_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_ne_en_s$bin_centers,
    density_diff_list_ne_en_s$density_diff 
  )

  # Twisted
  hist(EE_T,
       main="Twisted - EE",
       col=BAR_COLOR,
       border=twisted_color, 
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_ee_t)
  
  hist(NE_T,
       main="Twisted - NE (left)",
       col=BAR_COLOR,
       border=twisted_color, 
       xlab="Interaction distance",
       xlim=c(-XMAX,0), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topleft", cex=0.8, bty="n",
         legend=legend_vec_ne_t)
  
  hist(EN_T,
       main="Twisted - EN (right)",
       col=BAR_COLOR,
       border=twisted_color, 
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_en_t)
  
  hist(NN_T,
       main="Twisted - NN",
       col=BAR_COLOR,
       border=twisted_color,
       xlab="Interaction distance",
       xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
  legend("topright", cex=0.8, bty="n",
         legend=legend_vec_nn_t)
  
  # Plot differences of density differences of NE and EN within twisted
  plot(density_diff_list_ne_en_t$bin_centers,
       density_diff_list_ne_en_t$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - (NE-EN)",
       xlab="Interaction distance",
       ylab="NE minus EN",
       pch=21,
       col=BAR_COLOR,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_ne_en_t$bin_centers,
    density_diff_list_ne_en_t$density_diff 
  )
  
  # Get density differences of DI simple and DI twisted
  density_diff_list_ee_di_s_t <- function.get_density_diff(EE_S,EE_T, BREAKS)
  density_diff_list_ne_di_s_t <- function.get_density_diff(NE_S,NE_T, BREAKS)
  density_diff_list_en_di_s_t <- function.get_density_diff(EN_S,EN_T, BREAKS)
  density_diff_list_nn_di_s_t <- function.get_density_diff(NN_S,NN_T, BREAKS)
  
  MINY <- min(density_diff_list_ee_di_s_t$density_diff,
              density_diff_list_ne_di_s_t$density_diff,
              density_diff_list_en_di_s_t$density_diff,
              density_diff_list_nn_di_s_t$density_diff
  )
  
  MAXY <- max(density_diff_list_ee_di_s_t$density_diff,
              density_diff_list_ne_di_s_t$density_diff,
              density_diff_list_en_di_s_t$density_diff,
              density_diff_list_nn_di_s_t$density_diff
  )
  
  plot(density_diff_list_ee_di_s_t$bin_centers,
       density_diff_list_ee_di_s_t$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - EE",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_ee_di_s_t$bin_centers,
    density_diff_list_ee_di_s_t$density_diff 
  )
  
  plot(density_diff_list_ne_di_s_t$bin_centers,
       density_diff_list_ne_di_s_t$density_diff,
       xlim=c(-XMAX,0),
       ylim=c(MINY,MAXY),
       main="Density difference - NE",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=-270600, col="gray", lwd=0.3)
  abline(v=-270600*2, col="gray", lwd=0.3)
  abline(v=-270600*3, col="gray", lwd=0.3)
  abline(v=-270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_ne_di_s_t$bin_centers,
    density_diff_list_ne_di_s_t$density_diff 
  )
  
  plot(density_diff_list_en_di_s_t$bin_centers,
       density_diff_list_en_di_s_t$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - EN",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_en_di_s_t$bin_centers,
    density_diff_list_en_di_s_t$density_diff 
  )
  
  plot(density_diff_list_nn_di_s_t$bin_centers,
       density_diff_list_nn_di_s_t$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - NN",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_nn_di_s_t$bin_centers,
    density_diff_list_nn_di_s_t$density_diff 
  )
  
  plot.new()
  plot.new()
  
  # Plot dennsity different for NE and EN again with other ylim
  MINY <- min(density_diff_list_ne_di_s_t$density_diff,
              density_diff_list_en_di_s_t$density_diff
  )
  
  MAXY <- max(density_diff_list_ne_di_s_t$density_diff,
              density_diff_list_en_di_s_t$density_diff
  )
  
  plot(density_diff_list_ne_di_s_t$bin_centers,
       density_diff_list_ne_di_s_t$density_diff,
       xlim=c(-XMAX,0),
       ylim=c(MINY,MAXY),
       main="Density difference - NE",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=-270600, col="gray", lwd=0.3)
  abline(v=-270600*2, col="gray", lwd=0.3)
  abline(v=-270600*3, col="gray", lwd=0.3)
  abline(v=-270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_ne_di_s_t$bin_centers,
    density_diff_list_ne_di_s_t$density_diff 
  )
  
  plot(density_diff_list_en_di_s_t$bin_centers,
       density_diff_list_en_di_s_t$density_diff,
       xlim=c(0,XMAX),
       ylim=c(MINY,MAXY),
       main="Density difference - EN",
       xlab="Interaction distance",
       ylab="Simple minus twisted",
       pch=21,
       col=simple_color,
       bg=twisted_color)
  abline(h=0, col="gray")
  abline(v=270600, col="gray", lwd=0.3)
  abline(v=270600*2, col="gray", lwd=0.3)
  abline(v=270600*3, col="gray", lwd=0.3)
  abline(v=270600*4, col="gray", lwd=0.3)
  function.add_spline_curve(
    density_diff_list_en_di_s_t$bin_centers,
    density_diff_list_en_di_s_t$density_diff 
  )
  
  plot.new()
  
  # Add main title
  mtext(MAIN_MAIN_TITLE, outer = TRUE, cex = 1)
  
  dev.off()
  
}

function.write_statistics_to_file <- function(
  TSV_NAME,
  OUT_PREFIX,
  DI_EE_S,
  DI_NE_S,
  DI_EN_S,
  DI_NN_S,
  DI_EE_T,
  DI_NE_T,
  DI_EN_T,
  DI_NN_T,
  DI_EE_S_RP,
  DI_NE_S_RP,
  DI_EN_S_RP,
  DI_NN_S_RP,
  DI_EE_T_RP,
  DI_NE_T_RP,
  DI_EN_T_RP,
  DI_NN_T_RP,
  UIR_EE_S_RP,
  UIR_NE_S_RP,
  UIR_EN_S_RP,
  UIR_NN_S_RP,
  UIR_EE_T_RP,
  UIR_NE_T_RP,
  UIR_EN_T_RP,
  UIR_NN_T_RP,
  UI_EE_S_RP,
  UI_NE_S_RP,
  UI_EN_S_RP,
  UI_NN_S_RP,
  UI_EE_T_RP,
  UI_NE_T_RP,
  UI_EN_T_RP,
  UI_NN_T_RP
)
{
  # Calculate statistics
  # --------------------
  
  # Interaction numbers
  
  DI_EE_S_N <- length(DI_EE_S)
  DI_NE_S_N <- length(DI_NE_S)
  DI_EN_S_N <- length(DI_EN_S)
  DI_NN_S_N <- length(DI_NN_S)
  
  DI_EE_T_N <- length(DI_EE_T)
  DI_NE_T_N <- length(DI_NE_T)
  DI_EN_T_N <- length(DI_EN_T)
  DI_NN_T_N <- length(DI_NN_T)
  
  # Read pair numbers  
  
  DI_EE_S_RP_N <- length(DI_EE_S_RP)
  DI_NE_S_RP_N <- length(DI_NE_S_RP)
  DI_EN_S_RP_N <- length(DI_EN_S_RP)
  DI_NN_S_RP_N <- length(DI_NN_S_RP)
  
  DI_EE_T_RP_N <- length(DI_EE_T_RP)
  DI_NE_T_RP_N <- length(DI_NE_T_RP)
  DI_EN_T_RP_N <- length(DI_EN_T_RP)
  DI_NN_T_RP_N <- length(DI_NN_T_RP)
  
  UIR_EE_S_RP_N <- length(UIR_EE_S_RP)
  UIR_NE_S_RP_N <- length(UIR_NE_S_RP)
  UIR_EN_S_RP_N <- length(UIR_EN_S_RP)
  UIR_NN_S_RP_N <- length(UIR_NN_S_RP)
  
  UIR_EE_T_RP_N <- length(UIR_EE_T_RP)
  UIR_NE_T_RP_N <- length(UIR_NE_T_RP)
  UIR_EN_T_RP_N <- length(UIR_EN_T_RP)
  UIR_NN_T_RP_N <- length(UIR_NN_T_RP)
  
  UI_EE_S_RP_N <- length(UI_EE_S_RP)
  UI_NE_S_RP_N <- length(UI_NE_S_RP)
  UI_EN_S_RP_N <- length(UI_EN_S_RP)
  UI_NN_S_RP_N <- length(UI_NN_S_RP)
  
  UI_EE_T_RP_N <- length(UI_EE_T_RP)
  UI_NE_T_RP_N <- length(UI_NE_T_RP)
  UI_EN_T_RP_N <- length(UI_EN_T_RP)
  UI_NN_T_RP_N <- length(UI_NN_T_RP)
  
  # Median interaction distances
  
  DI_EE_S_MED <- round(median(DI_EE_S))
  DI_NE_S_MED <- round(median(DI_NE_S))
  DI_EN_S_MED <- round(median(DI_EN_S))
  DI_NN_S_MED <- round(median(DI_NN_S))
  
  DI_EE_T_MED <- round(median(DI_EE_T))
  DI_NE_T_MED <- round(median(DI_NE_T))
  DI_EN_T_MED <- round(median(DI_EN_T))
  DI_NN_T_MED <- round(median(DI_NN_T))
  
  # Median read pair distances  
  
  DI_EE_S_RP_MED <- round(median(DI_EE_S_RP))
  DI_NE_S_RP_MED <- round(median(DI_NE_S_RP))
  DI_EN_S_RP_MED <- round(median(DI_EN_S_RP))
  DI_NN_S_RP_MED <- round(median(DI_NN_S_RP))
  
  DI_EE_T_RP_MED <- round(median(DI_EE_T_RP))
  DI_NE_T_RP_MED <- round(median(DI_NE_T_RP))
  DI_EN_T_RP_MED <- round(median(DI_EN_T_RP))
  DI_NN_T_RP_MED <- round(median(DI_NN_T_RP))
  
  UIR_EE_S_RP_MED <- round(median(UIR_EE_S_RP))
  UIR_NE_S_RP_MED <- round(median(UIR_NE_S_RP))
  UIR_EN_S_RP_MED <- round(median(UIR_EN_S_RP))
  UIR_NN_S_RP_MED <- round(median(UIR_NN_S_RP))
  
  UIR_EE_T_RP_MED <- round(median(UIR_EE_T_RP))
  UIR_NE_T_RP_MED <- round(median(UIR_NE_T_RP))
  UIR_EN_T_RP_MED <- round(median(UIR_EN_T_RP))
  UIR_NN_T_RP_MED <- round(median(UIR_NN_T_RP))
  
  UI_EE_S_RP_MED <- round(median(UI_EE_S_RP))
  UI_NE_S_RP_MED <- round(median(UI_NE_S_RP))
  UI_EN_S_RP_MED <- round(median(UI_EN_S_RP))
  UI_NN_S_RP_MED <- round(median(UI_NN_S_RP))
  
  UI_EE_T_RP_MED <- round(median(UI_EE_T_RP))
  UI_NE_T_RP_MED <- round(median(UI_NE_T_RP))
  UI_EN_T_RP_MED <- round(median(UI_EN_T_RP))
  UI_NN_T_RP_MED <- round(median(UI_NN_T_RP))
  
  # Interquartile range of interaction distances
  
  DI_EE_S_IQR <- round(IQR(DI_EE_S))
  DI_NE_S_IQR <- round(IQR(DI_NE_S))
  DI_EN_S_IQR <- round(IQR(DI_EN_S))
  DI_NN_S_IQR <- round(IQR(DI_NN_S))
  
  DI_EE_T_IQR <- round(IQR(DI_EE_T))
  DI_NE_T_IQR <- round(IQR(DI_NE_T))
  DI_EN_T_IQR <- round(IQR(DI_EN_T))
  DI_NN_T_IQR <- round(IQR(DI_NN_T))
  
  # Interquartile range of read pair distances  
  
  DI_EE_S_RP_IQR <- round(IQR(DI_EE_S_RP))
  DI_NE_S_RP_IQR <- round(IQR(DI_NE_S_RP))
  DI_EN_S_RP_IQR <- round(IQR(DI_EN_S_RP))
  DI_NN_S_RP_IQR <- round(IQR(DI_NN_S_RP))
  
  DI_EE_T_RP_IQR <- round(IQR(DI_EE_T_RP))
  DI_NE_T_RP_IQR <- round(IQR(DI_NE_T_RP))
  DI_EN_T_RP_IQR <- round(IQR(DI_EN_T_RP))
  DI_NN_T_RP_IQR <- round(IQR(DI_NN_T_RP))
  
  UIR_EE_S_RP_IQR <- round(IQR(UIR_EE_S_RP))
  UIR_NE_S_RP_IQR <- round(IQR(UIR_NE_S_RP))
  UIR_EN_S_RP_IQR <- round(IQR(UIR_EN_S_RP))
  UIR_NN_S_RP_IQR <- round(IQR(UIR_NN_S_RP))
  
  UIR_EE_T_RP_IQR <- round(IQR(UIR_EE_T_RP))
  UIR_NE_T_RP_IQR <- round(IQR(UIR_NE_T_RP))
  UIR_EN_T_RP_IQR <- round(IQR(UIR_EN_T_RP))
  UIR_NN_T_RP_IQR <- round(IQR(UIR_NN_T_RP))
  
  UI_EE_S_RP_IQR <- round(IQR(UI_EE_S_RP))
  UI_NE_S_RP_IQR <- round(IQR(UI_NE_S_RP))
  UI_EN_S_RP_IQR <- round(IQR(UI_EN_S_RP))
  UI_NN_S_RP_IQR <- round(IQR(UI_NN_S_RP))
  
  UI_EE_T_RP_IQR <- round(IQR(UI_EE_T_RP))
  UI_NE_T_RP_IQR <- round(IQR(UI_NE_T_RP))
  UI_EN_T_RP_IQR <- round(IQR(UI_EN_T_RP))
  UI_NN_T_RP_IQR <- round(IQR(UI_NN_T_RP))

    
  # Write table line with all statistics to text file
  # -------------------------------------------------
  
  # Header line
  cat(c(
    
    "OUT_PREFIX",
    
    # Interaction numbers
    
    "DI_EE_S_N",
    "DI_NE_S_N",
    "DI_EN_S_N",
    "DI_NN_S_N",
    
    "DI_EE_T_N",
    "DI_NE_T_N",
    "DI_EN_T_N",
    "DI_NN_T_N",
    
    # Read pair numbers
    
    "DI_EE_S_RP_N",
    "DI_NE_S_RP_N",
    "DI_EN_S_RP_N",
    "DI_NN_S_RP_N",
    
    "DI_EE_T_RP_N",
    "DI_NE_T_RP_N",
    "DI_EN_T_RP_N",
    "DI_NN_T_RP_N",
    
    "UIR_EE_S_RP_N",
    "UIR_NE_S_RP_N",
    "UIR_EN_S_RP_N",
    "UIR_NN_S_RP_N",
    
    "UIR_EE_T_RP_N",
    "UIR_NE_T_RP_N",
    "UIR_EN_T_RP_N",
    "UIR_NN_T_RP_N",
    
    "UI_EE_S_RP_N",
    "UI_NE_S_RP_N",
    "UI_EN_S_RP_N",
    "UI_NN_S_RP_N",
    
    "UI_EE_T_RP_N",
    "UI_NE_T_RP_N",
    "UI_EN_T_RP_N",
    "UI_NN_T_RP_N",

    # Median interaction distances
    
    "DI_EE_S_MED",
    "DI_NE_S_MED",
    "DI_EN_S_MED",
    "DI_NN_S_MED",
    
    "DI_EE_T_MED",
    "DI_NE_T_MED",
    "DI_EN_T_MED",
    "DI_NN_T_MED",
    
    # Median read pair distances
    
    "DI_EE_S_RP_MED",
    "DI_NE_S_RP_MED",
    "DI_EN_S_RP_MED",
    "DI_NN_S_RP_MED",
    
    "DI_EE_T_RP_MED",
    "DI_NE_T_RP_MED",
    "DI_EN_T_RP_MED",
    "DI_NN_T_RP_MED",
    
    "UIR_EE_S_RP_MED",
    "UIR_NE_S_RP_MED",
    "UIR_EN_S_RP_MED",
    "UIR_NN_S_RP_MED",
    
    "UIR_EE_T_RP_MED",
    "UIR_NE_T_RP_MED",
    "UIR_EN_T_RP_MED",
    "UIR_NN_T_RP_MED",
    
    "UI_EE_S_RP_MED",
    "UI_NE_S_RP_MED",
    "UI_EN_S_RP_MED",
    "UI_NN_S_RP_MED",
    
    "UI_EE_T_RP_MED",
    "UI_NE_T_RP_MED",
    "UI_EN_T_RP_MED",
    "UI_NN_T_RP_MED",
    
    # Interquartile range of interaction distances
    
    "DI_EE_S_IQR",
    "DI_NE_S_IQR",
    "DI_EN_S_IQR",
    "DI_NN_S_IQR",
    
    "DI_EE_T_IQR",
    "DI_NE_T_IQR",
    "DI_EN_T_IQR",
    "DI_NN_T_IQR",
    
    # Interquartile range of read pair distances
    
    "DI_EE_S_RP_IQR",
    "DI_NE_S_RP_IQR",
    "DI_EN_S_RP_IQR",
    "DI_NN_S_RP_IQR",
    
    "DI_EE_T_RP_IQR",
    "DI_NE_T_RP_IQR",
    "DI_EN_T_RP_IQR",
    "DI_NN_T_RP_IQR",
    
    "UIR_EE_S_RP_IQR",
    "UIR_NE_S_RP_IQR",
    "UIR_EN_S_RP_IQR",
    "UIR_NN_S_RP_IQR",
    
    "UIR_EE_T_RP_IQR",
    "UIR_NE_T_RP_IQR",
    "UIR_EN_T_RP_IQR",
    "UIR_NN_T_RP_IQR",
    
    "UI_EE_S_RP_IQR",
    "UI_NE_S_RP_IQR",
    "UI_EN_S_RP_IQR",
    "UI_NN_S_RP_IQR",
    
    "UI_EE_T_RP_IQR",
    "UI_NE_T_RP_IQR",
    "UI_EN_T_RP_IQR",
    "UI_NN_T_RP_IQR"
    
  ), file=TSV_NAME, sep="\t")
  cat("\n", file=TSV_NAME, append=TRUE)
  cat(c(
    
    OUT_PREFIX,
    
    # Interaction numbers
    
    DI_EE_S_N,
    DI_NE_S_N,
    DI_EN_S_N,
    DI_NN_S_N,
    
    DI_EE_T_N,
    DI_NE_T_N,
    DI_EN_T_N,
    DI_NN_T_N,
    
    # Read pair numbers
    
    DI_EE_S_RP_N,
    DI_NE_S_RP_N,
    DI_EN_S_RP_N,
    DI_NN_S_RP_N,
    
    DI_EE_T_RP_N,
    DI_NE_T_RP_N,
    DI_EN_T_RP_N,
    DI_NN_T_RP_N,
    
    UIR_EE_S_RP_N,
    UIR_NE_S_RP_N,
    UIR_EN_S_RP_N,
    UIR_NN_S_RP_N,
    
    UIR_EE_T_RP_N,
    UIR_NE_T_RP_N,
    UIR_EN_T_RP_N,
    UIR_NN_T_RP_N,
    
    UI_EE_S_RP_N,
    UI_NE_S_RP_N,
    UI_EN_S_RP_N,
    UI_NN_S_RP_N,
    
    UI_EE_T_RP_N,
    UI_NE_T_RP_N,
    UI_EN_T_RP_N,
    UI_NN_T_RP_N,
    
    # Median interaction distances
    
    DI_EE_S_MED,
    DI_NE_S_MED,
    DI_EN_S_MED,
    DI_NN_S_MED,
    
    DI_EE_T_MED,
    DI_NE_T_MED,
    DI_EN_T_MED,
    DI_NN_T_MED,
    
    # Median read pair distances
    
    DI_EE_S_RP_MED,
    DI_NE_S_RP_MED,
    DI_EN_S_RP_MED,
    DI_NN_S_RP_MED,
    
    DI_EE_T_RP_MED,
    DI_NE_T_RP_MED,
    DI_EN_T_RP_MED,
    DI_NN_T_RP_MED,
    
    UIR_EE_S_RP_MED,
    UIR_NE_S_RP_MED,
    UIR_EN_S_RP_MED,
    UIR_NN_S_RP_MED,
    
    UIR_EE_T_RP_MED,
    UIR_NE_T_RP_MED,
    UIR_EN_T_RP_MED,
    UIR_NN_T_RP_MED,
    
    UI_EE_S_RP_MED,
    UI_NE_S_RP_MED,
    UI_EN_S_RP_MED,
    UI_NN_S_RP_MED,
    
    UI_EE_T_RP_MED,
    UI_NE_T_RP_MED,
    UI_EN_T_RP_MED,
    UI_NN_T_RP_MED,
    
    # Interquartile range of interaction distances
    
    DI_EE_S_IQR,
    DI_NE_S_IQR,
    DI_EN_S_IQR,
    DI_NN_S_IQR,
    
    DI_EE_T_IQR,
    DI_NE_T_IQR,
    DI_EN_T_IQR,
    DI_NN_T_IQR,
    
    # Interquartile range of read pair distances
    
    DI_EE_S_RP_IQR,
    DI_NE_S_RP_IQR,
    DI_EN_S_RP_IQR,
    DI_NN_S_RP_IQR,
    
    DI_EE_T_RP_IQR,
    DI_NE_T_RP_IQR,
    DI_EN_T_RP_IQR,
    DI_NN_T_RP_IQR,
    
    UIR_EE_S_RP_IQR,
    UIR_NE_S_RP_IQR,
    UIR_EN_S_RP_IQR,
    UIR_NN_S_RP_IQR,
    
    UIR_EE_T_RP_IQR,
    UIR_NE_T_RP_IQR,
    UIR_EN_T_RP_IQR,
    UIR_NN_T_RP_IQR,
    
    UI_EE_S_RP_IQR,
    UI_NE_S_RP_IQR,
    UI_EN_S_RP_IQR,
    UI_NN_S_RP_IQR,
    
    UI_EE_T_RP_IQR,
    UI_NE_T_RP_IQR,
    UI_EN_T_RP_IQR,
    UI_NN_T_RP_IQR
    
  ), file=TSV_NAME, append=TRUE, sep="\t")
  cat("\n", file=TSV_NAME, append=TRUE)
}


PDF_NAME <- paste(OUT_DIR, OUT_PREFIX,"_di_distance_histograms_st_ee_ne_en_nn.pdf",sep="")
MM_TITLE <- paste(OUT_PREFIX," - Distances of simple and twisted interactions",sep="")
PDF_W <- 12.5
PDF_H <- 11
BAR_COLOR <- directed_color
function.plot_histograms_i_ee_ne_en_nn(
  PDF_NAME,
  PDF_W,
  PDF_H,
  MM_TITLE,
  BREAKS,
  XMAX,
  BAR_COLOR,
  DI_EE_S_DIST,
  DI_NE_S_DIST,
  DI_EN_S_DIST,
  DI_NN_S_DIST,
  DI_EE_T_DIST,
  DI_NE_T_DIST,
  DI_EN_T_DIST,
  DI_NN_T_DIST
)

PDF_NAME <- paste(OUT_DIR, OUT_PREFIX,"_di_rp_distance_histograms_st_ee_ne_en_nn.pdf",sep="")
MM_TITLE <- paste(OUT_PREFIX," - Distances of simple and twisted read pairs in directed interactions",sep="")
PDF_W <- 12.5
PDF_H <- 11
BAR_COLOR <- directed_color
function.plot_histograms_i_ee_ne_en_nn(
  PDF_NAME,
  PDF_W,
  PDF_H,
  MM_TITLE,
  BREAKS,
  XMAX,
  BAR_COLOR,
  DI_EE_S_RP_DIST,
  DI_NE_S_RP_DIST,
  DI_EN_S_RP_DIST,
  DI_NN_S_RP_DIST,
  DI_EE_T_RP_DIST,
  DI_NE_T_RP_DIST,
  DI_EN_T_RP_DIST,
  DI_NN_T_RP_DIST
)

PDF_NAME <- paste(OUT_DIR, OUT_PREFIX,"_uir_rp_distance_histograms_st_ee_ne_en_nn.pdf",sep="")
MM_TITLE <- paste(OUT_PREFIX," - Distances of simple and twisted read pairs in undirected reference interactions",sep="")
PDF_W <- 12.5
PDF_H <- 11
BAR_COLOR <- undirected_ref_color
function.plot_histograms_i_ee_ne_en_nn(
  PDF_NAME,
  PDF_W,
  PDF_H,
  MM_TITLE,
  BREAKS,
  XMAX,
  BAR_COLOR,
  UIR_EE_S_RP_DIST,
  UIR_NE_S_RP_DIST,
  UIR_EN_S_RP_DIST,
  UIR_NN_S_RP_DIST,
  UIR_EE_T_RP_DIST,
  UIR_NE_T_RP_DIST,
  UIR_EN_T_RP_DIST,
  UIR_NN_T_RP_DIST
)

PDF_NAME <- paste(OUT_DIR, OUT_PREFIX,"_ui_rp_distance_histograms_st_ee_ne_en_nn.pdf",sep="")
MM_TITLE <- paste(OUT_PREFIX," - Distances of simple and twisted read pairs in undirected interactions",sep="")
PDF_W <- 12.5
PDF_H <- 11
BAR_COLOR <- undirected_color
function.plot_histograms_i_ee_ne_en_nn(
  PDF_NAME,
  PDF_W,
  PDF_H,
  MM_TITLE,
  BREAKS,
  XMAX,
  BAR_COLOR,
  UI_EE_S_RP_DIST,
  UI_NE_S_RP_DIST,
  UI_EN_S_RP_DIST,
  UI_NN_S_RP_DIST,
  UI_EE_T_RP_DIST,
  UI_NE_T_RP_DIST,
  UI_EN_T_RP_DIST,
  UI_NN_T_RP_DIST
)

TSV_NAME <- paste(OUT_DIR, OUT_PREFIX,"_st_distance_statistics_ee_ne_en_nn.tsv",sep="")
function.write_statistics_to_file(
  TSV_NAME,
  OUT_PREFIX,
  DI_EE_S_DIST,
  DI_NE_S_DIST,
  DI_EN_S_DIST,
  DI_NN_S_DIST,
  DI_EE_T_DIST,
  DI_NE_T_DIST,
  DI_EN_T_DIST,
  DI_NN_T_DIST,
  DI_EE_S_RP_DIST,
  DI_NE_S_RP_DIST,
  DI_EN_S_RP_DIST,
  DI_NN_S_RP_DIST,
  DI_EE_T_RP_DIST,
  DI_NE_T_RP_DIST,
  DI_EN_T_RP_DIST,
  DI_NN_T_RP_DIST,
  UIR_EE_S_RP_DIST,
  UIR_NE_S_RP_DIST,
  UIR_EN_S_RP_DIST,
  UIR_NN_S_RP_DIST,
  UIR_EE_T_RP_DIST,
  UIR_NE_T_RP_DIST,
  UIR_EN_T_RP_DIST,
  UIR_NN_T_RP_DIST,
  UI_EE_S_RP_DIST,
  UI_NE_S_RP_DIST,
  UI_EN_S_RP_DIST,
  UI_NN_S_RP_DIST,
  UI_EE_T_RP_DIST,
  UI_NE_T_RP_DIST,
  UI_EN_T_RP_DIST,
  UI_NN_T_RP_DIST
)
