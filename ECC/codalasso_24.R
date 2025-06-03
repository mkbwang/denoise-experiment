
rm(list=ls())
library(BenchmarkDenoise)
library(phyloseq)


args <- commandArgs(trailingOnly=TRUE)
myseed <- as.integer(args[1])


phyasv_24 <- readRDS("ECC/phyasv_visit24.rds")
metadata_24 <- sample_data(phyasv_24)

target <- metadata_24$CaseEver
target_binary <- 1*(target == "Case")

count_24_normalized <- read.csv("ECC/denoised/count_24_marker_original.csv",
                                row.names=1) |> as.matrix()
count_24_normalized_imputed <- simple_impute(count_mat=count_24_normalized)



aucs_df <- data.frame(train_auc=rep(0, 7),
                      test_auc=rep(0, 7))

coefs_mat <- matrix(0, nrow=7, ncol=ncol(count_24_normalized_imputed))
rownames(aucs_df) <- rownames(coefs_mat) <- c("Original", "Denoise_4", "Denoise_5",
                         "Denoise_6", "Denoise_7", "Denoise_8", "Denoise_9")
colnames(coefs_mat) <- colnames(count_24_normalized_imputed)

# lasso on original data
lasso_result <- fit_codalasso(y=target_binary, X=count_24_normalized_imputed,
                              lambdas=0.4,
                              seed=myseed)

aucs_df$train_auc[1] <- lasso_result$train_auc
aucs_df$test_auc[1] <- lasso_result$test_auc
coefs_mat[1, ] <- lasso_result$coefs

# lasso on denoised data
for (j in 1:6){
  print(j)
  selected_id <- j+3
  denoised_result <- readRDS(sprintf("denoised/count_24_marker_denoised_%d.rds",
                                     selected_id))
  
  count_24_normalized_denoised <- denoised_result$denoised_counts
  rownames(count_24_normalized_denoised) <- rownames(count_24_normalized_imputed)
  colnames(count_24_normalized_denoised) <- colnames(count_24_normalized_imputed)
  
  lasso_result_denoised <- fit_codalasso(y=target_binary, X=count_24_normalized_denoised,
                                         lambdas=c(0.3, 0.4, 0.5),
                                         seed=myseed)
  
  aucs_df$train_auc[j+1] <- lasso_result_denoised$train_auc
  aucs_df$test_auc[j+1] <- lasso_result_denoised$test_auc
  coefs_mat[j+1, ] <- lasso_result_denoised$coefs
  
}

aucs_df$seed <- myseed
output <- list(aucs_df=aucs_df,
               coefs_mat=coefs_mat)

filename <- sprintf("performance_24_%d.rds", myseed)
saveRDS(output, file.path("codalasso", filename))


