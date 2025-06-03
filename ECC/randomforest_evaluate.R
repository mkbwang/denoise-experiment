

# first evaluate 24 month data

rm(list=ls())
library(BenchmarkDenoise)
library(phyloseq)
library(ggplot2)

phyasv_24 <- readRDS("ECC/phyasv_visit24.rds")
DAA_result_24 <- read.csv("ECC/DAA/DAA_result_yr2.csv", row.names=1)
marker_ASVs <- readRDS("ECC/DAA/marker_ASVs.rds")


metadata_24 <- sample_data(phyasv_24)

target <- metadata_24$CaseEver
target_binary <- 1*(target == "Case")

aucs_24_list <- list()
vimp_24_list <- list()
for (j in 1:100){
  
  result <- readRDS(sprintf("ECC/randomForest/performance_24_%d.rds", j))
  aucs_24_list[[j]] <- result$aucs_df
  vimp_24_list[[j]] <- result$vimp_mat
  
}

train_aucs <- lapply(aucs_24_list, function(x) x$train_auc)
train_aucs_matrix <- do.call(cbind, train_aucs)

test_aucs <- lapply(aucs_24_list, function(x) x$test_auc)
test_aucs_matrix <- do.call(cbind, test_aucs)


test_aucs_df <- t(test_aucs_matrix) |> as.data.frame()
colnames(test_aucs_df) <- rownames(aucs_24_list[[1]])

ggplot(test_aucs_df, aes(x=Original, y=Denoise_5)) +
  geom_point(size=1, alpha=0.7) + 
  geom_abline(intercept=0, slope=1, linetype="dashed", color="blue")+
  scale_y_continuous(limits=c(0.5, 1),
                     breaks=seq(0.5,1, 0.1))+
  scale_y_continuous(limits=c(0.5, 1),
                     breaks=seq(0.5, 1, 0.1))+
  xlab("Original Test AUC") + ylab("Denoised Test AUC")


wilcox.test(test_aucs_matrix[1,], test_aucs_matrix[7, ], paired=T)

mean(test_aucs_matrix[1,] < test_aucs_matrix[6, ])

View(coefs_24_list[[6]])




