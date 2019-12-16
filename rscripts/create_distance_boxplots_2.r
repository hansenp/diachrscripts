#x<-read.table("JAV_ERY_RALT_0.0018_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_ERY",".pdf", sep="")

#x<-read.table("JAV_MK_RALT_0.0019_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_MK",".pdf", sep="")

#x<-read.table("JAV_ACD4_RALT_0.0019_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_ACD4",".pdf", sep="")

#x<-read.table("JAV_EP_RALT_0.0017_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_EP",".pdf", sep="")

#x<-read.table("JAV_MAC_M0_RALT_0.0019_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_MAC_M0",".pdf", sep="")

#x<-read.table("JAV_MAC_M1_RALT_0.0019_interactions_with_genesymbols.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_MAC_M1",".pdf", sep="")

#x<-read.table("JAV_MK_RALT_0.0019_interactions_with_genesymbols.dc.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
#PDF_NAME <- paste("distance_boxplots_MK_DC",".pdf", sep="")

x<-read.table("JAV_ERY_RALT_0.0018_interactions_with_genesymbols.dc.tsv", na.strings="", colClasses=c("NULL","integer", "character", "NULL", "NULL", "character", "NULL"))
PDF_NAME <- paste("distance_boxplots_ERY_DC",".pdf", sep="")



cat("------------\n")

INDEFINABLE_ALL <- x[which(x[,2]=="NA"),1]
UNDIRECTED_ALL_ALL <- x[which(x[,2]=="U"),1]
UNDIRECTED_REF_ALL <- x[which(x[,2]=="URAA" | x[,2]=="URAI" | x[,2]=="URII"),1]
SIMPLE_ALL <- x[which(x[,2]=="S"),1]
TWISTED_ALL <- x[which(x[,2]=="T"),1]
DIRECTED_ALL <- c(SIMPLE_ALL,TWISTED_ALL)

N_INDEFINABLE_ALL <- length(INDEFINABLE_ALL)
N_UNDIRECTED_ALL_ALL <- length(UNDIRECTED_ALL_ALL)
N_UNDIRECTED_REF_ALL <- length(UNDIRECTED_REF_ALL)
N_SIMPLE_ALL <- length(SIMPLE_ALL)
N_TWISTED_ALL <- length(TWISTED_ALL)
N_DIRECTED_ALL <- length(DIRECTED_ALL)

print(paste("# N_INDEFINABLE_ALL:",N_INDEFINABLE_ALL))
print(paste("# N_UNDIRECTED_ALL_ALL:",N_UNDIRECTED_ALL_ALL))
print(paste("# N_UNDIRECTED_REF_ALL:",N_UNDIRECTED_REF_ALL))
print(paste("# N_SIMPLE_ALL:",N_SIMPLE_ALL))
print(paste("# N_TWISTED_ALL:",N_TWISTED_ALL))
print(paste("# N_DIRECTED_ALL:",N_DIRECTED_ALL))

cat("------------\n")

ENRICHMENT_FLAG <- "II"
INDEFINABLE_II <- x[which(x[,2]=="NA" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_ALL_II <- x[which(x[,2]=="U" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_REF_II <- x[which((x[,2]=="URAA" | x[,2]=="URAI" | x[,2]=="URII") & x[,3]==ENRICHMENT_FLAG),1]
SIMPLE_II <- x[which(x[,2]=="S" & x[,3]==ENRICHMENT_FLAG),1]
TWISTED_II <- x[which(x[,2]=="T" & x[,3]==ENRICHMENT_FLAG),1]
DIRECTED_II <- c(SIMPLE_II,TWISTED_II)

N_INDEFINABLE_II <- length(INDEFINABLE_II)
N_UNDIRECTED_ALL_II <- length(UNDIRECTED_ALL_II)
N_UNDIRECTED_REF_II <- length(UNDIRECTED_REF_II)
N_SIMPLE_II <- length(SIMPLE_II)
N_TWISTED_II <- length(TWISTED_II)
N_DIRECTED_II <- length(DIRECTED_II)

print(paste("# N_INDEFINABLE_II:",N_INDEFINABLE_II))
print(paste("# N_UNDIRECTED_ALL_II:",N_UNDIRECTED_ALL_II))
print(paste("# N_UNDIRECTED_REF_II:",N_UNDIRECTED_REF_II))
print(paste("# N_SIMPLE_II:",N_SIMPLE_II))
print(paste("# N_TWISTED_II:",N_TWISTED_II))
print(paste("# N_DIRECTED_II:",N_DIRECTED_II))

cat("------------\n")

ENRICHMENT_FLAG <- "AI"
INDEFINABLE_AI <- x[which(x[,2]=="NA" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_ALL_AI <- x[which(x[,2]=="U" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_REF_AI <- x[which((x[,2]=="URAA" | x[,2]=="URAI" | x[,2]=="URII") & x[,3]==ENRICHMENT_FLAG),1]
SIMPLE_AI <- x[which(x[,2]=="S" & x[,3]==ENRICHMENT_FLAG),1]
TWISTED_AI <- x[which(x[,2]=="T" & x[,3]==ENRICHMENT_FLAG),1]
DIRECTED_AI <- c(SIMPLE_AI,TWISTED_AI)

N_INDEFINABLE_AI <- length(INDEFINABLE_AI)
N_UNDIRECTED_ALL_AI <- length(UNDIRECTED_ALL_AI)
N_UNDIRECTED_REF_AI <- length(UNDIRECTED_REF_AI)
N_SIMPLE_AI <- length(SIMPLE_AI)
N_TWISTED_AI <- length(TWISTED_AI)
N_DIRECTED_AI <- length(DIRECTED_AI)

print(paste("# N_INDEFINABLE_AI:",N_INDEFINABLE_AI))
print(paste("# N_UNDIRECTED_ALL_AI:",N_UNDIRECTED_ALL_AI))
print(paste("# N_UNDIRECTED_REF_AI:",N_UNDIRECTED_REF_AI))
print(paste("# N_SIMPLE_AI:",N_SIMPLE_AI))
print(paste("# N_TWISTED_AI:",N_TWISTED_AI))
print(paste("# N_DIRECTED_AI:",N_DIRECTED_AI))

cat("------------\n")

ENRICHMENT_FLAG <- "AA"
INDEFINABLE_AA <- x[which(x[,2]=="NA" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_ALL_AA <- x[which(x[,2]=="U" & x[,3]==ENRICHMENT_FLAG),1]
UNDIRECTED_REF_AA <- x[which((x[,2]=="URAA" | x[,2]=="URAI" | x[,2]=="URII") & x[,3]==ENRICHMENT_FLAG),1]
SIMPLE_AA <- x[which(x[,2]=="S" & x[,3]==ENRICHMENT_FLAG),1]
TWISTED_AA <- x[which(x[,2]=="T" & x[,3]==ENRICHMENT_FLAG),1]
DIRECTED_AA <- c(SIMPLE_AA,TWISTED_AA)

N_INDEFINABLE_AA <- length(INDEFINABLE_AA)
N_UNDIRECTED_ALL_AA <- length(UNDIRECTED_ALL_AA)
N_UNDIRECTED_REF_AA <- length(UNDIRECTED_REF_AA)
N_SIMPLE_AA <- length(SIMPLE_AA)
N_TWISTED_AA <- length(TWISTED_AA)
N_DIRECTED_AA <- length(DIRECTED_AA)

print(paste("# N_INDEFINABLE_AA:",N_INDEFINABLE_AA))
print(paste("# N_UNDIRECTED_ALL_AA:",N_UNDIRECTED_ALL_AA))
print(paste("# N_UNDIRECTED_REF_AA:",N_UNDIRECTED_REF_AA))
print(paste("# N_SIMPLE_AA:",N_SIMPLE_AA))
print(paste("# N_TWISTED_AA:",N_TWISTED_AA))
print(paste("# N_DIRECTED_AA:",N_DIRECTED_AA))

cat("------------\n")


cairo_pdf(PDF_NAME, height=21, width=3)
par(mfrow=c(7,1))
options(scipen = -1)

# plot barchart with interaction numbers in different categories
bp<-c(N_UNDIRECTED_ALL_ALL, N_UNDIRECTED_REF_ALL, N_DIRECTED_ALL)
names(bp) <- as.character(c("U", "UR", "D"))
barplot(bp, main="All interactions", col="gray", ylab="Interactions")
bp<-c(N_UNDIRECTED_ALL_II, N_UNDIRECTED_REF_II, N_DIRECTED_II)
names(bp) <- as.character(c("U", "UR", "D"))
barplot(bp, main="II interactions", col="blue", ylab="Interactions")
bp<-c(N_UNDIRECTED_ALL_AI, N_UNDIRECTED_REF_AI, N_DIRECTED_AI)
names(bp) <- as.character(c("U", "UR", "D"))
barplot(bp, main="AI interactions", col="purple", ylab="Interactions")
bp<-c(N_UNDIRECTED_ALL_AA, N_UNDIRECTED_REF_AA, N_DIRECTED_AA)
names(bp) <- as.character(c("U", "UR", "D"))
barplot(bp, main="AA interactions", col="red", ylab="Interactions")

# plot indefinable distances for 'ALL', 'II', 'AI' and 'AA'
data<-list(
  INDEFINABLE_ALL,
  INDEFINABLE_II,
  INDEFINABLE_AI,
  INDEFINABLE_AA
)

labels <- c(
  "ALL",
  "II",
  "AI", 
  "AA"
)

boxplot(
  data,
  las = 2,
  col=c(
    "gray",
    "blue",
    "purple",
    "red"
  ),
  outline=FALSE,
  main="Indefinalble interactions",
  ylab="Digest distance (nt)",
  names=labels
)

# plot undirected vs. directed
data<-list(
  UNDIRECTED_ALL_ALL,
  UNDIRECTED_REF_ALL,
  DIRECTED_ALL,
  UNDIRECTED_ALL_II,
  UNDIRECTED_REF_II,
  DIRECTED_II,
  UNDIRECTED_ALL_AI,
  UNDIRECTED_REF_AI,
  DIRECTED_AI,
  UNDIRECTED_ALL_AA,
  UNDIRECTED_REF_AA,
  DIRECTED_AA
)

labels <- c(
  "U - ALL",
  "UR - ALL",
  "D -ALL",
  "U - II",
  "UR - II",
  "D -II",
  "U - AI",
  "UR - AI",
  "D -AI",
  "U - AA",
  "UR - AA",
  "D -AA"
)

boxplot(
  data,
  las = 2,
  col=c(
    "gray",
    "gray",
    "gray",
    "blue",
    "blue",
    "blue",
    "purple",
    "purple",
    "purple",
    "red",
    "red",
    "red"
  ),
  outline=FALSE,
  main="Undirected vs. directed interactions",
  ylab="Digest distance (nt)",
  names=labels
)

# plot simple vs. twisted
data<-list(
  SIMPLE_ALL,
  TWISTED_ALL,
  SIMPLE_II,
  TWISTED_II,
  SIMPLE_AI,
  TWISTED_AI,
  SIMPLE_AA,
  TWISTED_AA
)

labels <- c(
  "S - ALL",
  "T - ALL",
  "S - II",
  "T - II",
  "S - AI",
  "T - AI",
  "S - AA",
  "T - AA"
)

boxplot(
  data,
  las = 2,
  col=c(
    "gray",
    "gray",
    "blue",
    "blue",
    "purple",
    "purple",
    "red",
    "red"
  ),
  outline=FALSE,
  main="Simple vs. twisted",
  ylab="Digest distance (nt)",
  names=labels
)

#y <- seq(0, YMAX, 5000000)
#axis(2, at=y, labels = formatC(y, big.mark = ",", format = "d"), las = 2)

#end_point = 1 + length(labels) + length(labels)-1

#text(seq(1, end_point*0.5,by=1), par("usr")[3]-3, 
  #   srt = 45, adj= 1, xpd = TRUE,
  #   labels = labels, cex=0.85)

dev.off()

#WCT_s<-wilcox.test(UNDIRECTED_REF, SIMPLE, alternative = "two.sided", log=T, paired=F, exact=T)
#print(paste("P-value simple: ", WCT_s$p.value))

#M_t<-wilcox.test(UNDIRECTED_REF, TWISTED, alternative = "two.sided", log=T, paired=F)
#print(paste("P-value twisted: ", M_t$p.value))

#M_d<-wilcox.test(UNDIRECTED_REF, DIRECTED, alternative = "two.sided", log=T, paired=F)
#print(paste("P-value directed: ", M_d$p.value))
