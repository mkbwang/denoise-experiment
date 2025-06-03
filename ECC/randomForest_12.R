
rm(list=ls())
library(BenchmarkDenoise)
library(phyloseq)


args <- commandArgs(trailingOnly=TRUE)
myseed <- as.integer(args[1])

phyasv_12 <- readRDS("phyasv_visit12.rds")
metadata_12 <- sample_data(phyasv_12)

target <- metadata_12$CaseEver
target_binary <- 1*(target == "Case")

count_12_normalized <- read.csv("denoised/count_12_marker_original.csv",
                                row.names=1) |> as.matrix()
count_12_normalized_imputed <- simple_impute(count_mat=count_12_normalized)



aucs_df <- data.frame(train_auc=rep(0, 7),
                      test_auc=rep(0, 7))

vimp_mat <- matrix(0, nrow=7, ncol=ncol(count_12_normalized_imputed))
rownames(aucs_df) <- rownames(vimp_mat) <- c("Original", "Denoise_4", "Denoise_5",
                         "Denoise_6", "Denoise_7", "Denoise_8", "Denoise_9")
colnames(vimp_mat) <- colnames(count_12_normalized_imputed)

# lasso on original data
rf_result <- fit_rf(y=target, X=count_12_normalized_imputed,
                              seed=myseed)

aucs_df$train_auc[1] <- rf_result$train_auc
aucs_df$test_auc[1] <- rf_result$test_auc
vimp_mat[1, ] <- rf_result$variable_importance

# lasso on denoised data
for (j in 1:6){
  print(j)
  selected_id <- j+3
  denoised_result <- readRDS(sprintf("denoised/count_12_marker_denoised_%d.rds",
                                     selected_id))
  
  count_12_normalized_denoised <- denoised_result$denoised_counts
  rownames(count_12_normalized_denoised) <- rownames(count_12_normalized_imputed)
  colnames(count_12_normalized_denoised) <- colnames(count_12_normalized_imputed)
  
  rf_result_denoised <- fit_rf(y=target, X=count_12_normalized_denoised,
                                         seed=myseed)
  
  aucs_df$train_auc[j+1] <- rf_result_denoised$train_auc
  aucs_df$test_auc[j+1] <- rf_result_denoised$test_auc
  vimp_mat[j+1, ] <- rf_result_denoised$variable_importance
  
}

aucs_df$seed <- myseed
output <- list(aucs_df=aucs_df,
               vimp_mat=vimp_mat)

filename <- sprintf("performance_12_%d.rds", myseed)
saveRDS(output, file.path("randomForest", filename))



