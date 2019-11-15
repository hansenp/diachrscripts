
# define function that returns the SSE
calcSSE <- function(x, y, d, s ){
  loessMod <- try(loess(y ~ x, data=d, span=s), silent=T)
  res <- try(loessMod$residuals, silent=T)
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 0)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}

get_fitted_line <- function(DATA,XMAX){
  x <- sort(DATA[which(DATA[,1]<XMAX),1])
  y <- DATA[order(x),2] 
  loess_fit <- loess(y~x, span=SPAN)
  y <- predict(loess_fit)
  return(list(x=x,y=y))
}

########################################################################################################################

MIFSUD=T
SCHOENEFELDER=F
CHESI_BMP2=F
CHESI_HEPG2=F

if(MIFSUD) {
  print("Loading data for Mifsud ...")
  SIMPLE<-read.table("../mifsud_alt_digest_distances_simple.txt")
  TWISTED<-read.table("../mifsud_alt_digest_distances_twisted.txt")
  DIRECTED<-rbind(SIMPLE,TWISTED)
  UNDIRECTED_REF<-read.table("../mifsud_alt_digest_distances_undirected_reference.txt")
  UNDIRECTED<-read.table("../mifsud_alt_digest_distances_undirected.txt")
  #INDEFINABLE<-read.table("../mifsud_alt_digest_distances_indefinable.txt")
  print("... done.")
}

if(SCHOENEFELDER) {
  print("Loading data for Schönefelder ...")
  SIMPLE<-read.table("../schoenefelder_alt_distances_digest_distances_simple.txt")
  TWISTED<-read.table("../schoenefelder_alt_distances_digest_distances_twisted.txt")
  DIRECTED<-rbind(SIMPLE,TWISTED)
  UNDIRECTED_REF<-read.table("../schoenefelder_alt_distances_digest_distances_undirected_reference.txt")
  # UNDIRECTED<-read.table("../schoenefelder_alt_distances_digest_distances_undirected.txt")
  # INDEFINABLE<-read.table("../schoenefelder_alt_distances_digest_distances_indefinable.txt")
  print("... done.")
}

if(CHESI_BMP2) {
  print("Loading data for Chesi-BMP2 ...")
  SIMPLE<-read.table("../chesi_bmp2_alt_digest_distances_simple.txt")
  TWISTED<-read.table("../chesi_bmp2_alt_digest_distances_twisted.txt")
  DIRECTED<-rbind(SIMPLE,TWISTED)
  UNDIRECTED_REF<-read.table("../chesi_bmp2_alt_digest_distances_undirected_reference.txt")
  # UNDIRECTED<-read.table("../chesi_bmp2_alt_digest_distances_undirected.txt")
  # INDEFINABLE<-read.table("../chesi_bmp2_alt_digest_distances_indefinable.txt")
  print("... done.")
}

if(CHESI_HEPG2) {
  print("Loading data for Chesi-HEPG2 ...")
  SIMPLE<-read.table("../chesi_hepg2_alt_digest_distances_simple.txt")
  TWISTED<-read.table("../chesi_hepg2_alt_digest_distances_twisted.txt")
  DIRECTED<-rbind(SIMPLE,TWISTED)
  UNDIRECTED_REF<-read.table("../chesi_hepg2_alt_digest_distances_undirected_reference.txt")
  # UNDIRECTED<-read.table("../chesi_hepg2_alt_digest_distances_undirected.txt")
  # INDEFINABLE<-read.table("../chesi_hepg2_alt_digest_distances_indefinable.txt")
  print("... done.")
}

cairo_pdf("foo.pdf", width=5, height=7.5)
par(mfrow=c(3,2), las=1)

XMIN <- 0
XMAX <- 6000000#max(SIMPLE[,1], TWISTED[,1], UNDIRECTED_REF[,1], UNDIRECTED[,1])
XLAB <- "Digest distance"

YMIN <- 0#min(SIMPLE[,2], TWISTED[,2], UNDIRECTED_REF[,2], UNDIRECTED[,2])
YMAX <- max(SIMPLE[,2], TWISTED[,2], UNDIRECTED_REF[,2], UNDIRECTED[,2])
YLAB <- "Number of read pairs"

SPAN <- 0.75
NBIN <- c(100,1500)
SHOW_OUTLIERS <- 50

print("Creating scatterplot for simple interactions ...")
smoothScatter(SIMPLE, main="Simple", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nbin=NBIN, nrpoints=SHOW_OUTLIERS)
L_s <- get_fitted_line(DATA=SIMPLE, XMAX=XMAX)
lines(L_s$x, L_s$y, col = "red")
print("... done.")

print("Creating scatterplot for twisted interactions ...")
smoothScatter(TWISTED, main="Twisted", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nbin=NBIN, nrpoints=SHOW_OUTLIERS)
L_t <- get_fitted_line(DATA=TWISTED, XMAX=XMAX)
lines(L_t$x, L_t$y, col = "red")
print("... done.")

print("Creating scatterplot for directed interactions ...")
smoothScatter(DIRECTED, main="Directed", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nbin=NBIN, nrpoints=SHOW_OUTLIERS)
L_d <- get_fitted_line(DATA=DIRECTED, XMAX=XMAX)
lines(L_d$x, L_d$y, col = "red")
print("... done.")

print("Creating scatterplot for undirected reference interactions ...")
smoothScatter(UNDIRECTED_REF, main="Undirected reference", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nbin=NBIN, nrpoints=SHOW_OUTLIERS)
L_ur <- get_fitted_line(DATA=UNDIRECTED_REF, XMAX=XMAX)
lines(L_ur$x, L_ur$y, col = "red")
print("... done.")

print("Creating plot with fitted lines for direcred and undirected reference interactions ...")
YMAX<-max(L_d$y,L_ur$y)
XMAX<-max(L_d$x,L_ur$x)
plot(L_d$x, L_d$y, col = "black", type="l", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab="Digest distance", ylab="Read pairs")
lines(L_ur$x, L_ur$y, col = "blue")
print("... done.")

#smoothScatter(UNDIRECTED, main="Undirected", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0, nbin=NBIN)
#L <- get_fitted_line(DATA=UNDIRECTED, XMAX=XMAX)
#lines(L$x, L$y, col = "red")

#plot(INDEFINABLE, main="Indefinable") # too large

dev.off()

