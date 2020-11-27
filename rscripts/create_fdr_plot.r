args = commandArgs(trailingOnly=TRUE)

P_THRESH=args[1]
PROTOCOL_COLOR=args[2]
PREFIX=args[3]
CSV_INPUT_FILE_LIST=args[4]

REPLICATE_COLORS <- c("blue","green","red","orange")
XMIN <- 0.000
XMAX <- 0.005
YMIN <- 0.000
YMAX <- 0.25
YMIN2 <- 0
YMAX2 <- 1000000
  
print(paste("P_THRESH:", P_THRESH))
print(paste("PROTOCOL_COLOR:", PROTOCOL_COLOR)) # one color for each protocol: "darkblue",capture Hi-C; "darkgreen",Capture-C; "chocolate4",Hi-C
print(paste("PREFIX:", PREFIX))
print(paste("CSV_INPUT_FILE_LIST:", CSV_INPUT_FILE_LIST)) # comma separted list of input files

# determine number of replicates
FILE_LIST <- unlist(strsplit(CSV_INPUT_FILE_LIST, ","))
REP_NUM <- length(FILE_LIST)-1
cat("\n")
print(paste("Number of replicates:",REP_NUM))

# read files to data frame
FDR <- vector()
PC <- vector()
NSIG <- vector()
for(f in FILE_LIST){
  print(paste("Reading file: ", f))
  x <- read.table(f, header=T)
  FDR <- cbind(FDR, x[,2])
  print(head(x[,2]))
  PC <- cbind(PC, x[,3])
  NSIG <- cbind(NSIG, x[,5])
}


cairo_pdf(paste(PREFIX, "_fdr_plot.pdf",sep=""), width=4.75, height=4)
par(mar = c(5,5,2,6), las=1)
for(i in 1:REP_NUM){
  plot(PC[,i], FDR[,i], xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), axes=F, xlab="", ylab="", col=REPLICATE_COLORS[i], pch=20)
  par(new=T)
}
plot(PC[,i+1], FDR[,i+1], xlim=c(XMIN,XMAX), ylim=c(YMIN,YMAX), xlab="P-value threshold", ylab="Estimated FDR", xaxp=c(0, 0.005, 2), yaxp=c(0, 1, 10), col=PROTOCOL_COLOR, pch=3)
abline(h=0.05, lty=2)
par(new=T)
plot(PC[,i+1], NSIG[,i+1], xlim=c(0,XMAX), ylim=c(0,YMAX2), axes=F, ann=FALSE, pch=18, col=PROTOCOL_COLOR)
axis(side = 4)
mtext(side = 4, line = 4, 'Selected interactions',las=0)
abline(v=P_THRESH, lty=2)
  
# Add a legend
PCH_VECTOR <- rep(20,REP_NUM)
PCH_VECTOR <- c(PCH_VECTOR,3,18)

print(REPLICATE_COLORS[1:REP_NUM])
print(PCH_VECTOR)

legend("topleft", legend=c(paste("Replicate", 1:REP_NUM, sep=" "),"At least two","At least two (selected)"), col=c(REPLICATE_COLORS[1:REP_NUM],PROTOCOL_COLOR,PROTOCOL_COLOR), pch=PCH_VECTOR, cex=0.8)
dev.off()

stop()

#REP_NUM=3

#PREFIX <- "../results_2/mifsud/mifsud"
#P_THRESH=0.00375
#PROTOCOL_COLOR="darkblue" # capture Hi-C

#PREFIX <- "../results_2/chesi/bmp2/chesi_bmp2"
#P_THRESH=0.0015
#PROTOCOL_COLOR="darkgreen" # Capture-C

#PREFIX <- "../results_2/chesi/hepg2/chesi_hepg2"
#P_THRESH=0.0015
#PROTOCOL_COLOR="darkgreen" # Capture-C

REP_NUM=2

PREFIX <- "../results_2/schoenefelder/schoenefelder"
P_THRESH=0.00175
PROTOCOL_COLOR="darkblue" # capture Hi-C


#PREFIX <- "../results_2/nora/untreated/nora_untreated"
#PREFIX <- "../results_2/nora/treated/nora_treated"
#PREFIX <- "../results_2/nora/washoff/nora_washoff"
#PREFIX <- "../results_2/schoenefelder/schoenefelder"



f_name<-paste(PREFIX, "_fdr_plot.pdf", sep="")
print(f_name)
cairo_pdf(f_name, width=4.75, height=4)

f_name<-paste(PREFIX, "_r1_fdr_analysis_results.txt", sep="")
M_R1 <- read.table(f_name, header=T)
m_r1_fdr <- M_R1[,2]
m_r1_pc <- M_R1[,3]
m_r1_nsigo <-M_R1[,5]

f_name<-paste(PREFIX, "_r2_fdr_analysis_results.txt", sep="")
M_R2 <- read.table(f_name, header=T)
m_r2_fdr <- M_R2[,2]
m_r2_pc <- M_R2[,3]
m_r2_nsigo <-M_R2[,5]

if(REP_NUM==3){
  f_name<-paste(PREFIX, "_r3_fdr_analysis_results.txt", sep="")
  M_R3 <- read.table(f_name, header=T)
  m_r3_fdr <- M_R3[,2]
  m_r3_pc <- M_R3[,3]
  m_r3_nsigo <-M_R3[,5]
}

f_name<-paste(PREFIX, "_alt_fdr_analysis_results.txt", sep="")
M_ALT <- read.table(f_name, header=T)
m_alt_fdr <- M_ALT[,2]
m_alt_pc <- M_ALT[,3]
m_alt_nsigo <-M_ALT[,5]


XMAX <- 0.005
YMAX <- 1.0
YMAX2 <- 100000

par(mar = c(5,5,2,6), las=1)
plot(m_r1_pc, m_r1_fdr, col="blue", ylim=c(0,YMAX), xlim=c(0,XMAX), xaxt='n', yaxt='n', ann=FALSE)
par(new=T)
plot(m_r2_pc, m_r2_fdr, col="green", ylim=c(0,YMAX), xlim=c(0,XMAX), xaxt='n', yaxt='n', ann=FALSE)
if(REP_NUM==3){
  par(new=T)
  plot(m_r3_pc, m_r3_fdr, col="red", xlim=c(0,XMAX), ylim=c(0,YMAX), xaxt='n', yaxt='n', ann=FALSE)
}
par(new=T)
plot(m_alt_pc, m_alt_fdr, col=PROTOCOL_COLOR, pch=3, xlim=c(0,XMAX), ylim=c(0,YMAX), xlab="P-value threshold", ylab="Estimated FDR", xaxp  = c(0, 0.005, 2), yaxp=c(0,1,4))
abline(h=0.25, lty=2)
par(new=T)
plot(m_alt_pc,m_alt_nsigo, xlim=c(0,XMAX), ylim=c(0,YMAX2), axes=F, ann=FALSE, pch=15, col=PROTOCOL_COLOR)
axis(side = 4)
mtext(side = 4, line = 4, 'Selected interactions')
abline(v=P_THRESH, lty=2)

# Add a legend
if(REP_NUM==3){
  legend("topleft", legend=c("Replicate 1", "Replicate 2", "Replicate 3", "At least two", "At least two (selected)"),
       col=c("blue", "green", "red", PROTOCOL_COLOR, PROTOCOL_COLOR), pch=c(1,1,1,3,15), cex=0.8)
} else {
  legend("topleft", legend=c("Replicate 1", "Replicate 2", "At least two", "At least two (selected)"),
         col=c("blue", "green", PROTOCOL_COLOR, PROTOCOL_COLOR), pch=c(1,1,3,15), cex=0.8)  
}
dev.off()