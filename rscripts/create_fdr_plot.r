REP_NUM=3
PREFIX <- "chesi_bmp2"
PREFIX <- "mifsud"
REP_NUM=2
PREFIX <- "nora_untreated"
PREFIX <- "schoenefelder"

f_name<-paste(PREFIX, "_fdr_plot.pdf", sep="")
cairo_pdf(f_name, width=6, height=5.5)

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
m_alt_nsigo <-M_ALT[,5]*m_alt_fdr

XMAX=max(m_r1_pc, m_r2_pc, m_r3_pc, m_alt_pc)
XMAX <- 0.005
YMAX=max(m_r1_fdr, m_r2_fdr, m_r3_fdr, m_alt_fdr)
YMAX <- 1.0

par(mar = c(5,5,2,5))
plot(m_r1_pc, m_r1_fdr, col="blue", ylim=c(0,YMAX), xlim=c(0,XMAX), xaxt='n', yaxt='n', ann=FALSE)
par(new=T)
plot(m_r2_pc, m_r2_fdr, col="green", ylim=c(0,YMAX), xlim=c(0,XMAX), xaxt='n', yaxt='n', ann=FALSE)
if(REP_NUM==3){
  par(new=T)
  plot(m_r3_pc, m_r3_fdr, col="red", xlim=c(0,XMAX), ylim=c(0,YMAX), xaxt='n', yaxt='n', ann=FALSE)
}
par(new=T)
plot(m_alt_pc, m_alt_fdr, col="black", pch=3, xlim=c(0,XMAX), ylim=c(0,YMAX), xlab="P-value cutoff", ylab="Estimated FDR")
abline(h=0.25, lty=2)
par(new=T)
plot(m_alt_pc,m_alt_nsigo, xlim=c(0,XMAX), ylim=c(0,25000), axes=F, ann=FALSE, pch=15)
axis(side = 4)
mtext(side = 4, line = 3, 'Number of selected interactions')

# Add a legend
if(REP_NUM==3){
  legend("topleft", legend=c("Replicate 1", "Replicate 2", "Replicate 3", "At least two", "At least two (selected)"),
       col=c("blue", "green", "red", "black", "black"), pch=c(1,1,1,3,15), cex=0.8)
} else {
  legend("topleft", legend=c("Replicate 1", "Replicate 2", "At least two", "At least two (selected)"),
         col=c("blue", "green", "black", "black"), pch=c(1,1,3,15), cex=0.8)  
}
dev.off()