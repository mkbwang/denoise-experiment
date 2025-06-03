
rm(list=ls())
library(DFBM)
library(BenchmarkDenoise)
library(circlize)


args <- commandArgs(trailingOnly=TRUE)
myseed <- as.integer(args[1])

myseed <- 20
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

col_fun <- colorRamp2(c(0, 1), c("#FFFFFF", "#222222"))

save_heatmap(X=mean_prob_mat, entry_name="Probability",
             filename="matrix_factorization/probability_mat.png",
             cmap=col_fun,
             width=20, height=6)

set.seed(myseed)
sim_vals <- matrix(rbinom(n=length(mean_prob_mat), size=1, prob=mean_prob_mat),
                   nrow=nrow(mean_prob_mat),
                   ncol=ncol(mean_prob_mat))

save_heatmap(X=sim_vals, entry_name="Probability",
             filename="matrix_factorization/simulated_mat.png",
             cmap=col_fun,
             width=20, height=6)


for (j in 1:7){
  
  sim_vals_missing <- sim_vals
  missing_ratio <- j/10
  missing_ID <- sample(length(sim_vals), length(sim_vals)*missing_ratio)
  sim_vals_missing[missing_ID] <- NA
  
  save_heatmap(X=sim_vals_missing, entry_name="Probability",
               filename=sprintf("matrix_factorization/simulated_mat_%d.png", j),
               cmap=col_fun,
               na_col="#8888F6",
               width=20, height=6)
  
}








