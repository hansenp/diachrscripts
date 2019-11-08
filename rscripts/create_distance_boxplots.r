
PREFIX <- "../results/mifsud/digest_distance/digest_distance_analysis_mifsud"

M_SIMPLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_simple.txt", sep=""))
M_TWISTED <- read.table(paste(PREFIX, "_at_least_2_digest_distances_twisted.txt", sep=""))
M_UNDIRECTED_REF <- read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected_reference.txt", sep=""))
M_UNDIRECTED_ALL <-read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected.txt", sep=""))
M_INDEFINABLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_indefinable.txt", sep=""))

PREFIX <- "../results/schoenefelder/digest_distance/digest_distance_analysis_schoenefelder"

S_SIMPLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_simple.txt", sep=""))
S_TWISTED <- read.table(paste(PREFIX, "_at_least_2_digest_distances_twisted.txt", sep=""))
S_UNDIRECTED_REF <- read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected_reference.txt", sep=""))
S_UNDIRECTED_ALL <-read.table(paste(PREFIX, "_at_least_2_digest_distances_undirected.txt", sep=""))
S_INDEFINABLE <- read.table(paste(PREFIX, "_at_least_2_digest_distances_indefinable.txt", sep=""))

PDF_NAME <- "_distance_plots_all.pdf"

print("Hurz")

data<-c(
  M_SIMPLE,
  M_TWISTED,
  M_UNDIRECTED_REF,
  M_UNDIRECTED_ALL,
  M_INDEFINABLE,
  S_SIMPLE,
  S_TWISTED,
  S_UNDIRECTED_REF,
  S_UNDIRECTED_ALL,
  S_INDEFINABLE  
)

labels <- c(
  "Simple (M)",
  "Twisted (M)",
  "Undirected ref. (M)",
  "Undirected all (M)",
  "Indefinable (M)",
  "Simple (S)",
  "Twisted (S)",
  "Undirected ref. (S)",
  "Undirected all (S)",
  "Indefinable (S)"
)

cairo_pdf("Hurz.pdf", width=6, height=5.5)
boxplot(
  data,
  las = 2,
  ylim=c(0,1000000),
  xaxt = "n",
  xlab = "",
  col=c(
    "orange",
    "orange",
    "orange",
    "orange",
    "orange",
    "lightblue",
    "lightblue",
    "lightblue",
    "lightblue",
    "lightblue"
  )
)


text(1:10, par("usr")[3] - 0.25, srt = 45, adj = 1, labels = labels, xpd = TRUE)


dev.off()