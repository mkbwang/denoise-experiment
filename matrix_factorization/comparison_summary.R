

library(ggplot2)
library(dplyr)

results_list <- list()

for (j in 1:100){
  result <- read.csv(sprintf("matrix_factorization/output/result_%d.csv", j))
  result$run <- j
  results_list[[j]] <- result
}


all_results_df <- do.call(rbind,results_list)
all_results_df$missingness <- sprintf("Left-out Proportion=%.1f", all_results_df$missingness)

write.csv(all_results_df, "matrix_factorization/simulation_results.csv",
          row.names = FALSE)

subset_results_df <- all_results_df %>% 
  filter(missingness != "Left-out Proportion=0.8")

pvals <- data.frame(Missing=seq(0.1, 0.7, 0.1),
                    AUC = 0,
                    Corr = 0,
                    MAE = 0)

for (j in 1:7){
  selected_results <- subset_results_df %>% 
    filter(missingness == sprintf("Left-out Proportion=%.1f", j/10))
  
   wilcoxtest_AUC <- wilcox.test(selected_results$lsvd_bestAUC, 
                                   selected_results$lcf_bestAUC, paired=TRUE)
   pvals$AUC[j] <- wilcoxtest_AUC$p.value
   
   wilcoxtest_corr <- wilcox.test(selected_results$lsvd_bestcor, 
                                 selected_results$lcf_bestcor, paired=TRUE)
   pvals$Corr[j] <- wilcoxtest_corr$p.value
   
   wilcoxtest_MAE <- wilcox.test(selected_results$lsvd_bestMAE, 
                                  selected_results$lcf_bestMAE, paired=TRUE)
   pvals$MAE[j] <- wilcoxtest_MAE$p.value
  
}

write.csv(all_results_df, "experiment/matrix_factorization/pvals.csv",
          row.names = FALSE)


ggplot(subset_results_df, aes(x=lsvd_bestAUC, y=lcf_bestAUC)) + 
  geom_point(size=0.7, alpha=0.7) + geom_abline(slope=1, intercept=0, linetype="dashed")+
  facet_wrap(vars(missingness), nrow=2)+
  xlab("LogisticSVD") + ylab("LogisticCF") + ggtitle("AUROC of Left-out Entries")


ggplot(subset_results_df, aes(x=lsvd_bestcor, y=lcf_bestcor)) + 
  geom_point(size=0.7, alpha=0.7) + geom_abline(slope=1, intercept=0, linetype="dashed")+
  facet_wrap(vars(missingness), nrow=2)+
  xlab("LogisticSVD") + ylab("LogisticCF") + ggtitle("Correlation with True Probabilities for Left-out Entries")


ggplot(subset_results_df, aes(x=lsvd_bestMAE, y=lcf_bestMAE)) + 
  geom_point(size=0.7, alpha=0.7) + geom_abline(slope=1, intercept=0, linetype="dashed")+
  facet_wrap(vars(missingness), nrow=2)+
  xlab("LogisticSVD") + ylab("LogisticCF") + ggtitle("Median Absolute Error with True Probabilities for Left-out Entries")


