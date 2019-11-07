

SIMPLE<-read.table("../schoenefelder_alt_distances_digest_distances_simple.txt")
TWISTED<-read.table("../schoenefelder_alt_distances_digest_distances_twisted.txt")
DIRECTED<-rbind(SIMPLE,TWISTED)
UNDIRECTED_REF<-read.table("../schoenefelder_alt_distances_digest_distances_undirected_reference.txt")
UNDIRECTED<-read.table("../schoenefelder_alt_distances_digest_distances_undirected.txt")
INDEFINABLE<-read.table("../schoenefelder_alt_distances_digest_distances_indefinable.txt")



cairo_pdf("foo.pdf", width=5, height=7.5)
par(mfrow=c(3,2), las=1)

XMIN <- 0
XMAX <- 10000000#max(SIMPLE[,1], TWISTED[,1], UNDIRECTED_REF[,1], UNDIRECTED[,1])
XLAB <- "Digest distance"

YMIN <- -100
YMAX <- max(SIMPLE[,2], TWISTED[,2], UNDIRECTED_REF[,2], UNDIRECTED[,2])
YLAB <- "Number of read pairs"

smoothScatter(SIMPLE, main="Simple", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0)
x <- sort(SIMPLE[,1])
y <- SIMPLE[order(SIMPLE[,1]),2] 
loess_fit_s <- loess(y~x)
lines(x, predict(loess_fit_s), col = "red")

smoothScatter(TWISTED, main="Twisted", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0)
x <- sort(TWISTED[,1])
y <- TWISTED[order(TWISTED[,1]),2] 
loess_fit_t <- loess(y~x)
lines(x, predict(loess_fit_t), col = "red")

smoothScatter(DIRECTED, main="Directed", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0)
x <- sort(DIRECTED[,1])
y <- DIRECTED[order(DIRECTED[,1]),2] 
loess_fit_d <- loess(y~x)
lines(x, predict(loess_fit_d), col = "red")

smoothScatter(UNDIRECTED_REF, main="Undirected reference", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0)
x <- sort(UNDIRECTED_REF[,1])
y <- UNDIRECTED_REF[order(UNDIRECTED_REF[,1]),2] 
loess_fit_ur <- loess(y~x)
lines(x, predict(loess_fit_ur), col = "red")

smoothScatter(UNDIRECTED, main="Undirected", xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab=XLAB, ylab=YLAB, nrpoints=0)


plot(INDEFINABLE, main="Indefinable")

dev.off()

