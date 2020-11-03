#!/usr/bin/env Rscript
source("rscripts/07_analyze_interaction_distances/interaction_distances_lib.r")

# Read commannd line arguments
args = commandArgs(trailingOnly=TRUE)

for(arg in args){
  print(arg)
}

OUT_DIR <- args[1]
OUT_PREFIX <- args[2]

DI_EE_DIST <- t(read.table(args[3]))
DI_NE_DIST <- t(-read.table(args[4]))
DI_EN_DIST <- t(read.table(args[5]))
DI_NN_DIST <- t(read.table(args[6]))

UIR_EE_DIST <- t(read.table(args[7]))
UIR_NE_DIST <- t(-read.table(args[8]))
UIR_EN_DIST <- t(read.table(args[9]))
UIR_NN_DIST <- t(read.table(args[10]))

UI_EE_DIST <- t(read.table(args[11]))
UI_NE_DIST <- t(-read.table(args[12]))
UI_EN_DIST <- t(read.table(args[13]))
UI_NN_DIST <- t(read.table(args[14]))


# Define auxiliary functions
# --------------------------

function.plot_histograms_ee_ne_en_nn <- function(
  PDF_FILE_NAME,
  PDF_WIDTH,
  PDF_HEIGHT,
  MAIN_MAIN_TITLE,
  BREAKS,
  XMAX,
  DI_EE,
  DI_NE,
  DI_EN,
  DI_NN,
  UIR_EE,
  UIR_NE,
  UIR_EN,
  UIR_NN,
  UI_EE,
  UI_NE,
  UI_EN,
  UI_NN  
  ) 
{
  # Get legend vectors with interaction numbers, median distances and median interquartile ranges

  legend_vec_di_ee <- function.get_legend_vector(DI_EE)
  legend_vec_di_ne <- function.get_legend_vector(DI_NE)
  legend_vec_di_en <- function.get_legend_vector(DI_EN)
  legend_vec_di_nn <- function.get_legend_vector(DI_NN)
  
  legend_vec_uir_ee <- function.get_legend_vector(UIR_EE)  
  legend_vec_uir_ne <- function.get_legend_vector(UIR_NE)
  legend_vec_uir_en <- function.get_legend_vector(UIR_EN)
  legend_vec_uir_nn <- function.get_legend_vector(UIR_NN)
  
  legend_vec_ui_ee <- function.get_legend_vector(UI_EE)  
  legend_vec_ui_ne <- function.get_legend_vector(UI_NE)
  legend_vec_ui_en <- function.get_legend_vector(UI_EN)
  legend_vec_ui_nn <- function.get_legend_vector(UI_NN)
  
  # Get common maximum for y-axes
  YMAX <- function.get_common_ymax_for_histograms(
    list(
      DI_EE,
      DI_NE,
      DI_EN,
      DI_NN,
      UIR_EE,
      UIR_NE,
      UIR_EN,
      UIR_NN,
      UI_EE,
      UI_NE,
      UI_EN,
      UI_NN      
      ),
    BREAKS)
  
  # Plot six histograms to PDF file
  cairo_pdf(PDF_FILE_NAME, width=PDF_WIDTH, height=PDF_HEIGHT)
  
    par(mfrow=c(9,5), oma = c(0, 0, 2, 0))

    # Directed
    hist(DI_EE,
         main="Directed - EE",
         col=directed_color,
         xlab="Interaction distance",
         xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
    legend("topright", cex=0.8, bty="n",
           legend=legend_vec_di_ee)
        
    hist(DI_NE,
         main="Directed - NE (left)",
         col=directed_color,
         xlab="Interaction distance",
         xlim=c(-XMAX,0), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topleft", cex=0.8, bty="n",
             legend=legend_vec_di_ne)
      
    hist(DI_EN,
         main="Directed - EN (right)",
         col=directed_color,
         xlab="Interaction distance",
         xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topright", cex=0.8, bty="n",
             legend=legend_vec_di_en)
      
      hist(DI_NN,
           main="Directed - NN",
           col=directed_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topright", cex=0.8, bty="n",
             legend=legend_vec_di_nn)
      
      # Get density differences of NE and EN for DI, UIR and UI
      density_diff_list_ee_dine_dien <- function.get_density_diff(-DI_NE,DI_EN, BREAKS)
      density_diff_list_ee_uirne_uiren <- function.get_density_diff(-UIR_NE,UIR_EN, BREAKS)
      density_diff_list_ee_uine_uien <- function.get_density_diff(-UI_NE,UI_EN, BREAKS)
      
      # Get common ylims of density differences for DI, UIR and UI
      MINY <- min(density_diff_list_ee_dine_dien$density_diff,
                  density_diff_list_ee_uirne_uiren$density_diff,
                  density_diff_list_ee_uine_uien$density_diff
      )
      MAXY <- max(density_diff_list_ee_dine_dien$density_diff,
                  density_diff_list_ee_uirne_uiren$density_diff,
                  density_diff_list_ee_uine_uien$density_diff
      )
      
      # Plot differences of density differences of NE and EN within DI
      plot(density_diff_list_ee_dine_dien$bin_centers,
           density_diff_list_ee_dine_dien$density_diff,
           col=directed_color,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - (NE-EN)",
           xlab="Interaction distance",
           ylab="NE minus EN",
           pch=20)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ee_dine_dien$bin_centers,
        density_diff_list_ee_dine_dien$density_diff
      )
      
      # Undirected reference
      hist(UIR_EE,
           main="Undirected reference - EE",
           col=undirected_ref_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topright", cex=0.8, bty="n",
             legend=legend_vec_uir_ee)
      
      hist(UIR_NE,
           main="Undirected reference - NE (left)",
           col=undirected_ref_color,
           xlab="Interaction distance",
           xlim=c(-XMAX,0), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
        legend("topleft", cex=0.8, bty="n",
               legend=legend_vec_uir_ne)
        
      hist(UIR_EN,
           main="Undirected reference - EN (right)",
           col=undirected_ref_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
        legend("topright", cex=0.8, bty="n",
               legend=legend_vec_uir_en)
      
      hist(UIR_NN,
           main="Undirected reference - NN",
           col=undirected_ref_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topright", cex=0.8, bty="n",
             legend=legend_vec_uir_nn)
      
      # Plot differences of density differences of NE and EN within UIR
      plot(density_diff_list_ee_uirne_uiren$bin_centers,
           density_diff_list_ee_uirne_uiren$density_diff,
           col=undirected_ref_color,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - (NE-EN)",
           xlab="Interaction distance",
           ylab="NE minus EN",
           pch=20)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ee_uirne_uiren$bin_centers,
        density_diff_list_ee_uirne_uiren$density_diff
      )

      # Undirected
      hist(UI_EE,
           main="Undirected - EE",
           col=undirected_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topright", cex=0.8, bty="n",
             legend=legend_vec_ui_ee)      
      
    hist(UI_NE,
         main="Undirected - NE (left)",
         col=undirected_color,
         xlab="Interaction distance",
         xlim=c(-XMAX,0), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
      legend("topleft", cex=0.8, bty="n",
             legend=legend_vec_ui_ne)
      
      hist(UI_EN,
           main="Undirected - EN (right)",
           col=undirected_color,
           xlab="Interaction distance",
           xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
        legend("topright", cex=0.8, bty="n",
               legend=legend_vec_ui_en)
  
        hist(UI_NN,
             main="Undirected - NN",
             col=undirected_color,
             xlab="Interaction distance",
             xlim=c(0,XMAX), breaks=BREAKS, ylim=c(0,YMAX), freq=F)
        legend("topright", cex=0.8, bty="n",
               legend=legend_vec_ui_nn)
        
        # Plot differences of density differences of NE and EN within UI
        plot(density_diff_list_ee_uine_uien$bin_centers,
             density_diff_list_ee_uine_uien$density_diff,
             col=undirected_color,
             xlim=c(0,XMAX),
             ylim=c(MINY,MAXY),
             main="Density difference - (NE-EN)",
             xlab="Interaction distance",
             ylab="NE minus EN",
             pch=20)
        abline(h=0, col="gray")
        abline(v=270600, col="gray", lwd=0.3)
        abline(v=270600*2, col="gray", lwd=0.3)
        abline(v=270600*3, col="gray", lwd=0.3)
        abline(v=270600*4, col="gray", lwd=0.3)
        function.add_spline_curve(
          density_diff_list_ee_uine_uien$bin_centers,
          density_diff_list_ee_uine_uien$density_diff
        )

      # Get density differences of DI and UIR
      density_diff_list_ee_di_uir <- function.get_density_diff(DI_EE,UIR_EE, BREAKS)
      density_diff_list_ne_di_uir <- function.get_density_diff(DI_NE,UIR_NE, BREAKS)
      density_diff_list_en_di_uir <- function.get_density_diff(DI_EN,UIR_EN, BREAKS)
      density_diff_list_nn_di_uir <- function.get_density_diff(DI_NN,UIR_NN, BREAKS)
      
      # Get density differences of DI and UI
      density_diff_list_ee_di_ui <- function.get_density_diff(DI_EE,UI_EE, BREAKS)
      density_diff_list_ne_di_ui <- function.get_density_diff(DI_NE,UI_NE, BREAKS)
      density_diff_list_en_di_ui <- function.get_density_diff(DI_EN,UI_EN, BREAKS)
      density_diff_list_nn_di_ui <- function.get_density_diff(DI_NN,UI_NN, BREAKS)
      
      # Get density differences of UIR and UI
      density_diff_list_ee_uir_ui <- function.get_density_diff(UIR_EE,UI_EE, BREAKS)
      density_diff_list_ne_uir_ui <- function.get_density_diff(UIR_NE,UI_NE, BREAKS)
      density_diff_list_en_uir_ui <- function.get_density_diff(UIR_EN,UI_EN, BREAKS)
      density_diff_list_nn_uir_ui <- function.get_density_diff(UIR_NN,UI_NN, BREAKS)
      
      MINY <- min(density_diff_list_ee_di_uir$density_diff,
                  density_diff_list_ne_di_uir$density_diff,
                  density_diff_list_en_di_uir$density_diff,
                  density_diff_list_nn_di_uir$density_diff,
                  density_diff_list_ee_di_ui$density_diff,
                  density_diff_list_ne_di_ui$density_diff,
                  density_diff_list_en_di_ui$density_diff,
                  density_diff_list_nn_di_ui$density_diff,
                  density_diff_list_ee_uir_ui$density_diff,
                  density_diff_list_ne_uir_ui$density_diff,
                  density_diff_list_en_uir_ui$density_diff,
                  density_diff_list_nn_uir_ui$density_diff
      )
      
      MAXY <- max(density_diff_list_ee_di_uir$density_diff,
                  density_diff_list_ne_di_uir$density_diff,
                  density_diff_list_en_di_uir$density_diff,
                  density_diff_list_nn_di_uir$density_diff,
                  density_diff_list_ee_di_ui$density_diff,
                  density_diff_list_ne_di_ui$density_diff,
                  density_diff_list_en_di_ui$density_diff,
                  density_diff_list_nn_di_ui$density_diff,
                  density_diff_list_ee_uir_ui$density_diff,
                  density_diff_list_ne_uir_ui$density_diff,
                  density_diff_list_en_uir_ui$density_diff,
                  density_diff_list_nn_uir_ui$density_diff
      )
      
      plot(density_diff_list_ee_di_uir$bin_centers,
           density_diff_list_ee_di_uir$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EE",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ee_di_uir$bin_centers,
        density_diff_list_ee_di_uir$density_diff
        )
      
      plot(density_diff_list_ne_di_uir$bin_centers,
           density_diff_list_ne_di_uir$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_di_uir$bin_centers,
        density_diff_list_ne_di_uir$density_diff
      )
      
      plot(density_diff_list_en_di_uir$bin_centers,
           density_diff_list_en_di_uir$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_di_uir$bin_centers,
        density_diff_list_en_di_uir$density_diff
      )
      
      
      
      plot(density_diff_list_nn_di_uir$bin_centers,
           density_diff_list_nn_di_uir$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - NN",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_nn_di_uir$bin_centers,
        density_diff_list_nn_di_uir$density_diff
      )
      
      plot.new()
      
      # Plot density differences of DI and UI
      plot(density_diff_list_ee_di_ui$bin_centers,
           density_diff_list_ee_di_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EE",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ee_di_ui$bin_centers,
        density_diff_list_ee_di_ui$density_diff
      )
      
      plot(density_diff_list_ne_di_ui$bin_centers,
           density_diff_list_ne_di_ui$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_di_ui$bin_centers,
        density_diff_list_ne_di_ui$density_diff
      )
      
      plot(density_diff_list_en_di_ui$bin_centers,
           density_diff_list_en_di_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_di_ui$bin_centers,
        density_diff_list_en_di_ui$density_diff
      )

      plot(density_diff_list_nn_di_ui$bin_centers,
           density_diff_list_nn_di_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - NN",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_nn_di_ui$bin_centers,
        density_diff_list_nn_di_ui$density_diff
      )
      
      plot.new()
      
      # Plot density differences of UIR and UI
      plot(density_diff_list_ee_uir_ui$bin_centers,
           density_diff_list_ee_uir_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EE",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ee_uir_ui$bin_centers,
        density_diff_list_ee_uir_ui$density_diff
      )
      
      plot(density_diff_list_ne_uir_ui$bin_centers,
           density_diff_list_ne_uir_ui$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_uir_ui$bin_centers,
        density_diff_list_ne_uir_ui$density_diff
      )
      
      plot(density_diff_list_en_uir_ui$bin_centers,
           density_diff_list_en_uir_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_uir_ui$bin_centers,
        density_diff_list_en_uir_ui$density_diff
      )

      plot(density_diff_list_nn_uir_ui$bin_centers,
           density_diff_list_nn_uir_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - NN",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_nn_uir_ui$bin_centers,
        density_diff_list_nn_uir_ui$density_diff
      )
      
      # Plot density differences again, but for NE and EN only
      MINY <- min(density_diff_list_ne_di_uir$density_diff,
                  density_diff_list_en_di_uir$density_diff,
                  density_diff_list_ne_di_ui$density_diff,
                  density_diff_list_en_di_ui$density_diff,
                  density_diff_list_ne_uir_ui$density_diff,
                  density_diff_list_en_uir_ui$density_diff
      )
      
      MAXY <- max(density_diff_list_ne_di_uir$density_diff,
                  density_diff_list_en_di_uir$density_diff,
                  density_diff_list_ne_di_ui$density_diff,
                  density_diff_list_en_di_ui$density_diff,
                  density_diff_list_ne_uir_ui$density_diff,
                  density_diff_list_en_uir_ui$density_diff
      )
      
      plot.new()
      plot.new()
      
      plot(density_diff_list_ne_di_uir$bin_centers,
           density_diff_list_ne_di_uir$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_di_uir$bin_centers,
        density_diff_list_ne_di_uir$density_diff
      )
      
      plot(density_diff_list_en_di_uir$bin_centers,
           density_diff_list_en_di_uir$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="DI minus UIR",
           pch=21,
           col=directed_color,
           bg=undirected_ref_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_di_uir$bin_centers,
        density_diff_list_en_di_uir$density_diff
      )
      
      plot.new()
      plot.new()
      
      plot.new()
      
      plot(density_diff_list_ne_di_ui$bin_centers,
           density_diff_list_ne_di_ui$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_di_ui$bin_centers,
        density_diff_list_ne_di_ui$density_diff
      )
      
      plot(density_diff_list_en_di_ui$bin_centers,
           density_diff_list_en_di_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="DI minus UI",
           pch=21,
           col=directed_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_di_ui$bin_centers,
        density_diff_list_en_di_ui$density_diff
      )
      
      plot.new()
      plot.new()
      
      plot.new()
      
      plot(density_diff_list_ne_uir_ui$bin_centers,
           density_diff_list_ne_uir_ui$density_diff,
           xlim=c(-XMAX,0),
           ylim=c(MINY,MAXY),
           main="Density difference - NE",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=-270600, col="gray", lwd=0.3)
      abline(v=-270600*2, col="gray", lwd=0.3)
      abline(v=-270600*3, col="gray", lwd=0.3)
      abline(v=-270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_ne_uir_ui$bin_centers,
        density_diff_list_ne_uir_ui$density_diff
      )
      
      plot(density_diff_list_en_uir_ui$bin_centers,
           density_diff_list_en_uir_ui$density_diff,
           xlim=c(0,XMAX),
           ylim=c(MINY,MAXY),
           main="Density difference - EN",
           xlab="Interaction distance",
           ylab="UIR minus UI",
           pch=21,
           col=undirected_ref_color,
           bg=undirected_color)
      abline(h=0, col="gray")
      abline(v=270600, col="gray", lwd=0.3)
      abline(v=270600*2, col="gray", lwd=0.3)
      abline(v=270600*3, col="gray", lwd=0.3)
      abline(v=270600*4, col="gray", lwd=0.3)
      function.add_spline_curve(
        density_diff_list_en_uir_ui$bin_centers,
        density_diff_list_en_uir_ui$density_diff
      )
      
      plot.new()
      plot.new()
      
    mtext(MAIN_MAIN_TITLE, outer = TRUE, cex = 1)
    
  dev.off()
}


# Write table line with all statistics to text file
function.write_statistics_to_file <- function(
  TSV_NAME,
  OUT_PREFIX,
  DI_EE,
  DI_NE,
  DI_EN,
  DI_NN,
  UIR_EE,
  UIR_NE,
  UIR_EN,
  UIR_NN,
  UI_EE,
  UI_NE,
  UI_EN,
  UI_NN 
)
{
  # Calculate statistics
  # --------------------
  
  # Interaction numbers
  DI_EE_N <- length(DI_EE)
  DI_NE_N <- length(DI_NE)
  DI_EN_N <- length(DI_EN)
  DI_NN_N <- length(DI_NN)
  UIR_EE_N <- length(UIR_EE)
  UIR_NE_N <- length(UIR_NE)
  UIR_EN_N <- length(UIR_EN)
  UIR_NN_N <- length(UIR_NN)
  UI_EE_N <- length(UI_EE)
  UI_NE_N <- length(UI_NE)
  UI_EN_N <- length(UI_EN)
  UI_NN_N <- length(UI_NN)
  
  # Median interaction distances
  DI_EE_MED <- round(median(DI_EE))
  DI_NE_MED <- round(median(DI_NE))
  DI_EN_MED <- round(median(DI_EN))
  DI_NN_MED <- round(median(DI_NN))
  UIR_EE_MED <- round(median(UIR_EE))
  UIR_NE_MED <- round(median(UIR_NE))
  UIR_EN_MED <- round(median(UIR_EN))
  UIR_NN_MED <- round(median(UIR_NN))
  UI_EE_MED <- round(median(UI_EE))
  UI_NE_MED <- round(median(UI_NE))
  UI_EN_MED <- round(median(UI_EN))
  UI_NN_MED <- round(median(UI_NN))
  
  # Interquartile ranges
  DI_EE_IQR <- round(IQR(DI_EE))
  DI_NE_IQR <- round(IQR(DI_NE))
  DI_EN_IQR <- round(IQR(DI_EN))
  DI_NN_IQR <- round(IQR(DI_NN))
  UIR_EE_IQR <- round(IQR(UIR_EE))
  UIR_NE_IQR <- round(IQR(UIR_NE))
  UIR_EN_IQR <- round(IQR(UIR_EN))
  UIR_NN_IQR <- round(IQR(UIR_NN))
  UI_EE_IQR <- round(IQR(UI_EE))
  UI_NE_IQR <- round(IQR(UI_NE))
  UI_EN_IQR <- round(IQR(UI_EN))
  UI_NN_IQR <- round(IQR(UI_NN))
  
  # Header line
  cat(c(
    
    "OUT_PREFIX",
    
    "DI_EE_N",
    "DI_NE_N",
    "DI_EN_N",
    "DI_NN_N",
    "UIR_EE_N",
    "UIR_NE_N",
    "UIR_EN_N",
    "UIR_NN_N",
    "UI_EE_N",
    "UI_NE_N",
    "UI_EN_N",
    "UI_NN_N",
    
    "DI_EE_MED",
    "DI_NE_MED",
    "DI_EN_MED",
    "DI_NN_MED",
    "UIR_EE_MED",
    "UIR_NE_MED",
    "UIR_EN_MED",
    "UIR_NN_MED",
    "UI_EE_MED",
    "UI_NE_MED",
    "UI_EN_MED",
    "UI_NN_MED",
    
    "DI_EE_IQR",
    "DI_NE_IQR",
    "DI_EN_IQR",
    "DI_NN_IQR",
    "UIR_EE_IQR",
    "UIR_NE_IQR",
    "UIR_EN_IQR",
    "UIR_NN_IQR",
    "UI_EE_IQR",
    "UI_NE_IQR",
    "UI_EN_IQR",
    "UI_NN_IQR"
    
  ), file=TSV_NAME, sep="\t")
  cat("\n", file=TSV_NAME, append=TRUE)
  cat(c(
    
    OUT_PREFIX,
    
    DI_EE_N,
    DI_NE_N,
    DI_EN_N,
    DI_NN_N,
    UIR_EE_N,
    UIR_NE_N,
    UIR_EN_N,
    UIR_NN_N,
    UI_EE_N,
    UI_NE_N,
    UI_EN_N,
    UI_NN_N,
    
    DI_EE_MED,
    DI_NE_MED,
    DI_EN_MED,
    DI_NN_MED,
    UIR_EE_MED,
    UIR_NE_MED,
    UIR_EN_MED,
    UIR_NN_MED,
    UI_EE_MED,
    UI_NE_MED,
    UI_EN_MED,
    UI_NN_MED,
    
    DI_EE_IQR,
    DI_NE_IQR,
    DI_EN_IQR,
    DI_NN_IQR,
    UIR_EE_IQR,
    UIR_NE_IQR,
    UIR_EN_IQR,
    UIR_NN_IQR,
    UI_EE_IQR,
    UI_NE_IQR,
    UI_EN_IQR,
    UI_NN_IQR
    
  ), file=TSV_NAME, append=TRUE, sep="\t")
  cat("\n", file=TSV_NAME, append=TRUE)
}


# Execute
# -------

# Plot histograms for all interaction categories
PDF_NAME <- paste(OUT_DIR, OUT_PREFIX,"_i_distance_histograms_ee_ne_en_nn.pdf",sep="")
MM_TITLE <- paste(OUT_PREFIX," - Interaction distances for DI, UIR and UI",sep="")
PDF_W <- 12.5
PDF_H <- 24.75
function.plot_histograms_ee_ne_en_nn(
  PDF_NAME,
  PDF_W,
  PDF_H,
  MM_TITLE,
  BREAKS,
  XMAX,
  DI_EE_DIST,
  DI_NE_DIST,
  DI_EN_DIST,
  DI_NN_DIST,
  UIR_EE_DIST,
  UIR_NE_DIST,
  UIR_EN_DIST,
  UIR_NN_DIST,
  UI_EE_DIST,
  UI_NE_DIST,
  UI_EN_DIST,
  UI_NN_DIST
)

# Write table with summary statistics (also shown in the legends) to a file for downstream analyzes
TSV_NAME <- paste(OUT_DIR, OUT_PREFIX,"_i_distance_statistics_ee_ne_en_nn.tsv",sep="")
function.write_statistics_to_file(
  TSV_NAME,
  OUT_PREFIX,
  DI_EE_DIST,
  -DI_NE_DIST,
  DI_EN_DIST,
  DI_NN_DIST,
  UIR_EE_DIST,
  -UIR_NE_DIST,
  UIR_EN_DIST,
  UIR_NN_DIST,
  UI_EE_DIST,
  -UI_NE_DIST,
  UI_EN_DIST,
  UI_NN_DIST 
)
