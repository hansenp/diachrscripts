
print("Reading distances for Mifsud ...")
PREFIX <- "../results_2/mifsud/mifsud_alt"

M_SIMPLE <- read.table(paste(PREFIX, "_digest_distances_simple.txt", sep=""))
M_SIMPLE <- M_SIMPLE[,1]
M_TWISTED <- read.table(paste(PREFIX, "_digest_distances_twisted.txt", sep=""))
M_TWISTED <- M_TWISTED[,1]
M_UNDIRECTED_REF <- read.table(paste(PREFIX, "_digest_distances_undirected_reference.txt", sep=""))
M_UNDIRECTED_REF <- M_UNDIRECTED_REF[,1]
M_UNDIRECTED_ALL <-read.table(paste(PREFIX, "_digest_distances_undirected.txt", sep=""))
M_UNDIRECTED_ALL <- M_UNDIRECTED_ALL[,1]
M_INDEFINABLE <- read.table(paste(PREFIX, "_digest_distances_indefinable.txt", sep=""))
M_INDEFINABLE <- M_INDEFINABLE[,1]


print("Reading distances for Schoenefelder ...")
PREFIX <- "../results_2/schoenefelder/schoenefelder_alt"

S_SIMPLE <- read.table(paste(PREFIX, "_digest_distances_simple.txt", sep=""))
S_SIMPLE <- S_SIMPLE[,1]
S_TWISTED <- read.table(paste(PREFIX, "_digest_distances_twisted.txt", sep=""))
S_TWISTED <- S_TWISTED[,1]
S_UNDIRECTED_REF <- read.table(paste(PREFIX, "_digest_distances_undirected_reference.txt", sep=""))
S_UNDIRECTED_REF <- S_UNDIRECTED_REF[,1]
S_UNDIRECTED_ALL <-read.table(paste(PREFIX, "_digest_distances_undirected.txt", sep=""))
S_UNDIRECTED_ALL <- S_UNDIRECTED_ALL[,1]
S_INDEFINABLE <- read.table(paste(PREFIX, "_digest_distances_indefinable.txt", sep=""))
S_INDEFINABLE <- S_INDEFINABLE[,1]

print("Creating boxplot ...")

PDF_NAME <- "distance_boxplots_all.pdf"

data<-list(
  M_INDEFINABLE,
  M_UNDIRECTED_ALL,
  M_UNDIRECTED_REF,
  M_SIMPLE,
  M_TWISTED,
  S_INDEFINABLE,
  S_UNDIRECTED_ALL,
  S_UNDIRECTED_REF,
  S_SIMPLE,
  S_TWISTED
)

labels <- c(
  "Indefinable (M)",
  "Undirected all (M)",
  "Undirected ref. (M)",
  "Simple (M)",
  "Twisted (M)",
  "Indefinable (S)",
  "Undirected all (S)",
  "Undirected ref. (S)",
  "Simple (S)",
  "Twisted (S)"
)

cairo_pdf(PDF_NAME, width=6, height=5.5)
#par(oma=c(0,0,0,0))
m<-quantile(M_INDEFINABLE)
s<-quantile(S_INDEFINABLE)

YMAX <- max(m[[4]],s[[4]])
boxplot(
  data,
  las = 2,
  ylim=c(0,YMAX),
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


#text(1:10, par("usr")[3] - 0.25, srt = 45, adj = 1, labels = labels, xpd = TRUE)
end_point = 1 + length(labels) + length(labels)-1

text(seq(1, end_point*0.5,by=1), par("usr")[3]-3, 
     srt = 45, adj= 1, xpd = TRUE,
     labels = labels, cex=1)


dev.off()