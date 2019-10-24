# excute from diachscripts/rscripts
PREFIX <- "../results_1/mifsud/mifsud_r1"
PREFIX <- "../results_1/mifsud/mifsud_r2"
PREFIX <- "../results_1/mifsud/mifsud_r3"
PREFIX <- "../results_1/mifsud/mifsud_alt"

PREFIX <- "../results_1/schoenefelder/schoenefelder_r1"
PREFIX <- "../results_1/schoenefelder/schoenefelder_r2"
PREFIX <- "../results_1/schoenefelder/schoenefelder_alt"

PREFIX <- "../results_1/chesi/bmp2/chesi_bmp2_r1"
PREFIX <- "../results_1/chesi/bmp2/chesi_bmp2_r2"
PREFIX <- "../results_1/chesi/bmp2/chesi_bmp2_r3"
PREFIX <- "../results_1/chesi/bmp2/chesi_bmp2_alt"

PREFIX <- "../results_1/chesi/hepg2/chesi_hepg2_r1"
PREFIX <- "../results_1/chesi/hepg2/chesi_hepg2_r2"
PREFIX <- "../results_1/chesi/hepg2/chesi_hepg2_r3"
PREFIX <- "../results_1/chesi/hepg2/chesi_hepg2_alt"

PREFIX <- "../results_1/nora/untreated/nora_untreated_r1"
PREFIX <- "../results_1/nora/untreated/nora_untreated_r2"
PREFIX <- "../results_1/nora/untreated/nora_untreated_alt"

PREFIX <- "../results_1/nora/treated/nora_treated_r1"
PREFIX <- "../results_1/nora/treated/nora_treated_r2"
#PREFIX <- "../results_1/nora/treated/nora_treated_alt"

#PREFIX <- "../results_1/nora/washoff/nora_washoff_r1"
#PREFIX <- "../results_1/nora/washoff/nora_washoff_r2"
#PREFIX <- "../results_1/nora/washoff/nora_washoff_alt"

PREFIX <- "../mifsud_alt_iter_1000"
#PREFIX <- "../schoenefelder_alt_iter_1000"
#PREFIX <- "../chesi_bmp2_alt_iter_1000"
PREFIX <- "../chesi_hepg2_alt_iter_1000"
PREFIX <- "../nora_untreated_alt_iter_1000"
PREFIX <- "../nora_treated_alt_iter_1000"
PREFIX <- "../nora_washoff_alt_iter_1000" 




f_name<-paste(PREFIX, "_permutation_plot.pdf", sep="")
cairo_pdf(f_name, width=4.5, height=4.5)

f_name<-paste(PREFIX, "_permutation_summary.txt", sep="")
SUMMARY<-read.table(f_name, header=T)

f_name<-paste(PREFIX, "_n_sig_permuted_interactions.txt", sep="")
INTERACTION_NUMBERS<-read.table(f_name)

OUT_PREFIX <- SUMMARY[,1]
INTERACTION_NUM <- SUMMARY[,2]
if(10000<=INTERACTION_NUM){
  INTERACTION_NUM <- format(as.numeric(INTERACTION_NUM), big.mark=",")
}
NOMINAL_ALPHA <- SUMMARY[,3]
ITER_NUM <- SUMMARY[,2]
NSIG_OBSERVED <- SUMMARY[,10]
NSIG_PERMUTED_MEAN <- round(SUMMARY[,11])
NSIG_RATIO_OBSERVED_PERMUTED <- round(NSIG_OBSERVED/NSIG_PERMUTED_MEAN,2)
BETTER_THAN_OBSERVED <- SUMMARY[,12]
N_TOTAL_INTERACION <- SUMMARY[,5]
N_INDEF_INTERACION <- SUMMARY[,8]
N_UNDIR_INTERACION <- SUMMARY[,9]

# compute Z score
p_mean <- mean(t(INTERACTION_NUMBERS))
p_sd <- sd(t(INTERACTION_NUMBERS))
z_score <- round((NSIG_OBSERVED-p_mean)/p_sd,2)
print(PREFIX)
print(paste("INTERACTION_NUM: ", INTERACTION_NUM))
print(paste("NSIG_OBSERVED: ", NSIG_OBSERVED))
print(paste("MEAN PERMUTED: ", round(p_mean,0)))
print(paste("SD PERMUTED: ", round(p_sd,2)))
print(paste("Z-score: ", z_score))

hist(t(INTERACTION_NUMBERS), xlim=c(min(INTERACTION_NUMBERS),NSIG_OBSERVED+20), main=OUT_PREFIX, xlab="Significant interactions")
abline(v=NSIG_OBSERVED, lty=2, col="red")

LEGEND_1 <- paste("Permuted interactions:", INTERACTION_NUM)
LEGEND_2 <- paste("Nominal alpha:", NOMINAL_ALPHA)
LEGEND_3 <- paste("Iterations:", ITER_NUM)
LEGEND_4 <- paste("Better than observerd:", BETTER_THAN_OBSERVED)
LEGEND_5 <- paste("Total number of interactions:", N_TOTAL_INTERACION)
LEGEND_6 <- paste("Number of indefinable interactions:", N_INDEF_INTERACION)
LEGEND_7 <- paste("Number of undirected interactions:", N_UNDIR_INTERACION)
LEGEND_8 <- paste("Number of significant directed interactions:", NSIG_OBSERVED)
LEGEND_9 <- paste("Mean permuted significant interactions:", NSIG_PERMUTED_MEAN)
LEGEND_10 <- paste("Ratio observed/permuted:", NSIG_RATIO_OBSERVED_PERMUTED)
LEGEND_11 <- paste("Z score:", z_score)
legend("topright", legend=c(LEGEND_1, LEGEND_2, LEGEND_3, LEGEND_4, LEGEND_5, LEGEND_6, LEGEND_7, LEGEND_8, LEGEND_9, LEGEND_10, LEGEND_11),
       col=c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"), pch=c(20,20,20,20,20,20,20,20,20,20,20), cex=0.6)

dev.off()
