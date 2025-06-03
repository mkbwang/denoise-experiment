
library(DFBM)
library(BenchmarkDenoise)
library(logisticPCA)
library(pROC)

expit <- function(x){
  1/(1+exp(-x))
}

args <- commandArgs(trailingOnly=TRUE)
myseed <- as.integer(args[1])

set.seed(myseed)

probs <- c(0.1, 0.3, 0.7, 0.9)
feature1 <- probs[sample(4, size=60, replace=T)]
feature2 <- probs[sample(4, size=60, replace=T)]
feature3 <- probs[sample(4, size=60, replace=T)]
feature4 <- probs[sample(4, size=60, replace=T)]

unique_features <- cbind(feature1, feature2, feature3, feature4)
qr(unique_features)$rank

sample_counts <- c(60, 60, 40, 40)

features1 <- matrix(rep(feature1, times=60), ncol=60, byrow=FALSE)
features2 <- matrix(rep(feature2, times=60), ncol=60, byrow=FALSE)
features3 <- matrix(rep(feature3, times=40), ncol=40, byrow=FALSE)
features4 <- matrix(rep(feature4, times=40), ncol=40, byrow=FALSE)


mean_prob_mat <- cbind(features1, features2, features3, features4)

set.seed(myseed)
sim_vals <- matrix(rbinom(n=length(mean_prob_mat), size=1, prob=mean_prob_mat),
                   nrow=nrow(mean_prob_mat),
                   ncol=ncol(mean_prob_mat))


comparison <- data.frame(missingness = seq(0.1, 0.8, 0.1),
                         lsvd_bestAUC = 0,
                         lsvd_bestcor = 0,
                         lsvd_bestMAE = 0,
                         lsvd_bestrank = 0,
                         lcf_bestAUC = 0,
                         lcf_bestcor = 0,
                         lcf_bestMAE = 0,
                         lcf_bestrank = 0)



# TODO: change seed
set.seed(myseed)

for (i in 1:8){
  
  missing_ratio <- comparison$missingness[i]
  
  validation_mask_ID <- sample(length(sim_vals), length(sim_vals)*missing_ratio)
  train_mask_ID <- setdiff(seq(1, length(sim_vals)), validation_mask_ID)
  
  binary_mask <- matrix(1, nrow=nrow(sim_vals), ncol=ncol(sim_vals))
  binary_mask[validation_mask_ID] <- 0
  
  # first apply logisticSVD
  sim_vals_withmissing <- sim_vals
  sim_vals_withmissing[validation_mask_ID] <- NA
  
  lsvd_prob_list <- list()
  for (selected_rank in seq(2,6)){
    lsvd_fit <- logisticSVD(sim_vals_withmissing, k=selected_rank-1, max_iters=5000)
    lsvd_logit <- t(replicate(nrow(sim_vals), lsvd_fit$mu))+
      lsvd_fit$A %*% t(lsvd_fit$B)
    lsvd_probs <- expit(lsvd_logit)
    lsvd_prob_list[[selected_rank-1]] <- lsvd_probs
  }
  
  lcf_prob_list <- list()
  for (selected_rank in seq(2,6)){
    # then try logisticCF
    lcf_fit <- logisticcfR(X=sim_vals, Z=binary_mask, K=selected_rank, lambda=0)
    lcf_logit <- t(lcf_fit$A) %*% lcf_fit$B
    lcf_probs <- expit(lcf_logit)
    lcf_prob_list[[selected_rank-1]] <- lcf_probs
  }
  
  
  AUC_train <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  AUC_validation <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  
  cor_train <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  cor_validation <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  
  mae_train <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  
  mae_validation <- data.frame(rank=seq(2,6), lsvd=rep(0, 5), lcf=rep(0, 5))
  
  
  for (j in 1:5){
    
    AUC_train$lsvd[j] <- auc(sim_vals[train_mask_ID], 
                             lsvd_prob_list[[j]][train_mask_ID])
    AUC_train$lcf[j] <- auc(sim_vals[train_mask_ID], 
                            lcf_prob_list[[j]][train_mask_ID])
    
    AUC_validation$lsvd[j] <- auc(sim_vals[validation_mask_ID], 
                                  lsvd_prob_list[[j]][validation_mask_ID])
    AUC_validation$lcf[j] <- auc(sim_vals[validation_mask_ID], 
                                 lcf_prob_list[[j]][validation_mask_ID])
    
    cor_train$lsvd[j] <- cor(lsvd_prob_list[[j]][train_mask_ID], 
                             mean_prob_mat[train_mask_ID])
    cor_train$lcf[j] <- cor(lcf_prob_list[[j]][train_mask_ID], 
                            mean_prob_mat[train_mask_ID])
    
    cor_validation$lsvd[j] <- cor(lsvd_prob_list[[j]][validation_mask_ID], 
                                  mean_prob_mat[validation_mask_ID])
    cor_validation$lcf[j] <- cor(lcf_prob_list[[j]][validation_mask_ID], 
                                 mean_prob_mat[validation_mask_ID])
    
    mae_train$lsvd[j] <- median(abs(lsvd_prob_list[[j]][train_mask_ID] - 
                                      mean_prob_mat[train_mask_ID]))
    mae_train$lcf[j] <- median(abs(lcf_prob_list[[j]][train_mask_ID] - 
                                     mean_prob_mat[train_mask_ID]))
    
    mae_validation$lsvd[j] <- median(abs(lsvd_prob_list[[j]][validation_mask_ID] - 
                                           mean_prob_mat[validation_mask_ID]))
    mae_validation$lcf[j] <- median(abs(lcf_prob_list[[j]][validation_mask_ID] - 
                                          mean_prob_mat[validation_mask_ID]))
    
  }
  
  comparison$lsvd_bestAUC[i] <- max(AUC_validation$lsvd)
  comparison$lsvd_bestcor[i] <- max(cor_validation$lsvd)
  comparison$lsvd_bestMAE[i] <- min(mae_validation$lsvd)
  comparison$lsvd_bestrank[i] <- AUC_validation$rank[which.max(AUC_validation$lsvd)]
  
  comparison$lcf_bestAUC[i] <- max(AUC_validation$lcf)
  comparison$lcf_bestcor[i] <- max(cor_validation$lcf)
  comparison$lcf_bestMAE[i] <- min(mae_validation$lcf)
  comparison$lcf_bestrank[i] <- AUC_validation$rank[which.max(AUC_validation$lcf)]
  
}


output_file <- sprintf("result_%d.csv", myseed)

write.csv(comparison, 
          file=file.path("output", output_file), 
          row.names=FALSE)



