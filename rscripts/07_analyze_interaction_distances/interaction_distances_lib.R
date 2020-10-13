# This file contains funnctions and variables that are used in more than one script

BIN_SIZE <- 20000
BREAKS <- seq(-300000000,300000000, BIN_SIZE)
XMAX <- 1000000

# Define colors
directed_color <- rgb(255/255,163/255,0/255,1)
undirected_ref_color <- rgb(171/255,215/255,230/255,1)
undirected_color <- rgb(210/255,210/255,210/255,1)
simple_color <- rgb(233/255,175/255,175/255,1)
twisted_color <- rgb(0/255,138/255,138/255,1)

function.get_legend_vector <- function(distance_vector)
{
  # Get interaction numbers
  n <- paste("n: ", format(length(distance_vector), big.mark=","), sep="")
  
  # Get median distances
  med <- paste("median: ", format(round(median(distance_vector)), big.mark=","), sep="")
  
  # Get median absolute deviation
  iqr <- paste("IQR: ", format(round(IQR(distance_vector)), big.mark=","), sep="")
  
  # Return legend vector
  return(c(n, med, iqr))
}

function.get_common_ymax_for_histograms <- function(distance_vectors_list, BREAKS)
{
  # Init ymax
  ym <- 0
  
  # Iterate list of distance vectors and update ymax
  for(d_vec in distance_vectors_list){
    d_vec <- hist(d_vec, breaks=BREAKS, plot=F)
    if(ym < max(d_vec$density)){ym <- max(d_vec$density)}
  }
  
  # Return ymax
  return(ym)
}

function.get_density_diff <- function(
  DIST_VEC_1,
  DIST_VEC_2,
  BREAKS
) 
{
  # Get densities in bins
  DIST_VEC_1 <- hist(DIST_VEC_1, breaks=BREAKS, plot=F)
  DIST_VEC_2 <- hist(DIST_VEC_2, breaks=BREAKS, plot=F)
  
  # Return bin centers shifted by the half of the bin size and density differeneces
  return(list(
    bin_centers=seq(BREAKS[1]+round((BREAKS[2]-BREAKS[1])/2),BREAKS[length(BREAKS)],(BREAKS[2]-BREAKS[1])),
    density_diff=DIST_VEC_1$density-DIST_VEC_2$density
  ))
}