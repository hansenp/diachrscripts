PDF_NAME <- "mifsud_distance_plots.pdf"
PREFIX <- "../results/mifsud/digest_distance/digest_distance_analysis_mifsud"

PDF_NAME <- "schoenefelder_distance_plots.pdf"
PREFIX <- "../results/schoenefelder/digest_distance/digest_distance_analysis_schoenefelder"


SIMPLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_simple.txt", sep=""))
TWISTED <- read.table(paste(PREFIX, "_at_least_2_digest_distances_twisted.txt", sep=""))
UNDIRECTED_REF <- read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected_reference.txt", sep=""))
UNDIRECTED_ALL <-read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected.txt", sep=""))
INDEFINABLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_indefinable.txt", sep=""))

DIRECTED <- c(t(SIMPLE), t(TWISTED))
print(nrow(UNDIRECTED_REF))
print(nrow(UNDIRECTED_ALL))
UNDIRECTED_ALL_SAMPLE <- sample(x=t(UNDIRECTED_ALL), size=nrow(UNDIRECTED_REF))


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




cairo_pdf(PDF_NAME, width=6, height=12)
par(mfrow=c(4,2))

XLAB="Distance (nt)"
YLIM <- c(0, max(H_SIMPLE$counts, H_TWISTED$counts))
hist(t(SIMPLE), breaks=BREAKS, xlim=XLIM, xlab=XLAB, main="Simple")
hist(t(TWISTED), breaks=BREAKS, xlim=XLIM, xlab=XLAB, main="Twisted")

YLIM <- c(0, max(H_DIRECTED$counts, H_UNDIRECTED_REF$counts))
hist(DIRECTED, breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Directed (simple and twisted)")
hist(t(UNDIRECTED_REF), breaks=BREAKS, xlim=XLIM, ylim=YLIM, xlab=XLAB, main="Undirected reference")


hist(t(UNDIRECTED_ALL), breaks=BREAKS, xlim=XLIM, xlab=XLAB, main="All undirected")
hist(t(INDEFINABLE), breaks=BREAKS,BIN_NUM, xlim=XLIM, xlab=XLAB, main="Indefinable")

YMIN <- min(H_DIRECTED$counts-H_UNDIRECTED_REF$counts, H_UNDIRECTED_REF$counts-H_UNDIRECTED_ALL_SAMPLE$counts)
YMAX <- max(H_DIRECTED$counts-H_UNDIRECTED_REF$counts, H_UNDIRECTED_REF$counts-H_UNDIRECTED_ALL_SAMPLE$counts)
YLIM=c(YMIN,YMAX)
plot(BREAKS, c(H_DIRECTED$counts-H_UNDIRECTED_REF$counts,0), xlim=c(0,1000000), ylim=YLIM, xlab=XLAB, ylab="Difference of bin counts", main="Directed minus undirected reference")
abline(h=0, lty=2)
plot(BREAKS, c(H_UNDIRECTED_REF$counts-H_UNDIRECTED_ALL_SAMPLE$counts,0), xlim=c(0,1000000), ylim=YLIM, xlab=XLAB, ylab="Difference of bin counts", main="Undirected ref. minus sampled undirected")
abline(h=0, lty=2)


dev.off()