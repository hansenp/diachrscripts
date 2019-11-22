PDF_NAME <- "mifsud_distance_plots_all.pdf"
PREFIX <- "../results_2/mifsud/mifsud_alt"
COLOR="orange"

PDF_NAME <- "schoenefelder_distance_plots_all.pdf"
PREFIX <- "../results_2/schoenefelder/schoenefelder_alt"
COLOR="lightblue"

#PDF_NAME <- "chesi_bmp2_distance_plots_all.pdf"
#PREFIX <- "../results_2/chesi/bmp2/chesi_bmp2_alt"
#COLOR="lightgreen"

#PDF_NAME <- "chesi_hepg2_distance_plots_all.pdf"
#PREFIX <- "../results_2/chesi/bmp2/chesi_hepg2_alt"
#COLOR="brown"

SIMPLE <- read.table(paste(PREFIX, "_digest_distances_simple.txt", sep=""))
SIMPLE <- SIMPLE[,1]
TWISTED <- read.table(paste(PREFIX, "_digest_distances_twisted.txt", sep=""))
TWISTED <- TWISTED[,1]
UNDIRECTED_REF <- read.table(paste(PREFIX, "_digest_distances_undirected_reference.txt", sep=""))
UNDIRECTED_REF <- UNDIRECTED_REF[,1]
UNDIRECTED_ALL <-read.table(paste(PREFIX, "_digest_distances_undirected.txt", sep=""))
UNDIRECTED_ALL <- UNDIRECTED_ALL[,1]
INDEFINABLE <- read.table(paste(PREFIX, "_digest_distances_indefinable.txt", sep=""))
INDEFINABLE <- INDEFINABLE[,1]

DIRECTED <- c(SIMPLE,TWISTED)
print(length(UNDIRECTED_REF))
print(length(UNDIRECTED_ALL))
UNDIRECTED_ALL_SAMPLE <- sample(x=t(UNDIRECTED_ALL), size=length(UNDIRECTED_REF))

LONGEST_INTERACTION <- max(SIMPLE,TWISTED,UNDIRECTED_REF,UNDIRECTED_ALL, INDEFINABLE)
XLIM <- c(0,1000000)

BIN_SIZE <- 20000
BREAKS <- seq(0, LONGEST_INTERACTION+BIN_SIZE, BIN_SIZE)

BIN_NUM <- LONGEST_INTERACTION/BIN_SIZE
print(LONGEST_INTERACTION/BIN_SIZE)

# determine YLIM of frequencies
H_SIMPLE <- hist(t(SIMPLE), breaks=BREAKS, plot=FALSE)
H_TWISTED <- hist(t(TWISTED), breaks=BREAKS, plot=FALSE)

H_DIRECTED <- hist(DIRECTED, breaks=BREAKS, plot=FALSE)
H_UNDIRECTED_REF <- hist(t(UNDIRECTED_REF), breaks=BREAKS, plot=FALSE)

H_UNDIRECTED_ALL <- hist(t(UNDIRECTED_ALL), breaks=BREAKS, plot=FALSE)
H_UNDIRECTED_ALL_SAMPLE <- hist(t(UNDIRECTED_ALL_SAMPLE), breaks=BREAKS, plot=FALSE)
H_INDEFINABLE <- hist(t(INDEFINABLE), breaks=BREAKS, plot=FALSE)

cairo_pdf(PDF_NAME, width=5, height=11)
par(mfrow=c(4,2))

XLAB="Distance (nt)"
YLIM <- c(0, max(H_SIMPLE$counts, H_TWISTED$counts))
hist(t(SIMPLE), breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Simple", col=COLOR)
hist(t(TWISTED), breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Twisted", col=COLOR)

YLIM <- c(0, max(H_DIRECTED$counts, H_UNDIRECTED_REF$counts))
hist(DIRECTED, breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Directed (simple and twisted)", col=COLOR)
hist(t(UNDIRECTED_REF), breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Undirected reference", col=COLOR)


hist(t(UNDIRECTED_ALL), breaks=BREAKS, xlim=XLIM, xlab=XLAB, main="All undirected", col=COLOR)
hist(t(INDEFINABLE), breaks=BREAKS,BIN_NUM, xlim=XLIM, xlab=XLAB, main="Indefinable", col=COLOR)

YMIN <- min(H_DIRECTED$counts-H_UNDIRECTED_REF$counts)
YMAX <- max(H_DIRECTED$counts-H_UNDIRECTED_REF$counts)
YLIM=c(YMIN,YMAX)
plot(BREAKS, c(H_DIRECTED$counts-H_UNDIRECTED_REF$counts,0), xlim=c(0,1000000), ylim=YLIM, xlab=XLAB, ylab="Difference between bin counts", main="Directed minus undirected reference", col=COLOR, pch=16, las=1)
abline(h=0, lty=2)
plot(BREAKS, c(H_UNDIRECTED_REF$counts-H_UNDIRECTED_ALL_SAMPLE$counts,0), xlim=c(0,1000000), ylim=YLIM, xlab=XLAB, ylab="Difference between bin counts", main="Undirected ref. minus sampled undirected", col=COLOR, pch=16, las=1)
abline(h=0, lty=2)

dev.off()