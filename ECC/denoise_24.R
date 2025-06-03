

rm(list=ls())

library(BenchmarkDenoise)
library(phyloseq)
library(caret)
library(ADAPT)
library(DFBM)

args <- commandArgs(trailingOnly=TRUE)
myseed <- as.integer(args[1])



# load phyloseq object
phyasv_24 <- readRDS("phyasv_visit24.rds")
metadata_24 <- sample_data(phyasv_24)
samples <- sample_names(phyasv_24)


# split into training and test samples
train_ids <- createDataPartition(metadata_24$CaseEver, p=0.6, list=FALSE)
train_samples <- samples[train_ids]
train_labels <- 1*(metadata_24$CaseEver[train_ids] == 'Case')
test_samples <- samples[-train_ids]
test_labels <- 1*(metadata_24$CaseEver[-train_ids] == 'Case')


# identify marker taxa
phyasv_24_train <- prune_samples(train_samples, phyasv_24)
count_train <- otu_table(phyasv_24_train)@.Data
phyasv_24_test <- prune_samples(test_samples, phyasv_24)
count_test <- otu_table(phyasv_24_test)@.Data
da_result_train <- adapt(phyasv_24_train, cond.var="CaseEver")
da_result_table_train <- summary(da_result_train)
marker_taxa <- da_result_table_train$Taxa[da_result_table_train$pval < 0.05]



# denoise training data and test data separately
count_train <- count_train[, marker_taxa]
count_test <- count_test[, marker_taxa]



count_train_normalized <- normalize(count_mat=count_train)
count_train_imputed <- simple_impute(count_mat=count_train_normalized) |> t()
count_test_normalized <- normalize(count_mat=count_test)
count_test_imputed <- simple_impute(count_mat=count_test_normalized) |> t()


output <- list()
output$train_labels <- train_labels
output$test_labels <- test_labels

output$train_original <- count_train_normalized
output$test_original <- count_test_normalized



for (j in 4:9){
  
  print(j)
  print("Denoise train")
  train_denoise <- dfbm(count_mat = count_train_normalized,
                        quantiles = seq(0.1, 0.9, 0.1),
                        increment=j/10,
                        max_K=6,
                        ignore=0.05,
                        interpolate=FALSE,
                        ncores=4)
  count_train_denoised <- train_denoise$denoised_counts
  rownames(count_train_denoised) <- rownames(count_train_normalized)
  colnames(count_train_denoised) <- colnames(count_train_normalized)
  output[[sprintf("train_denoised_%d", j)]] <- count_train_denoised
  
  
  print("Denoise test")
  test_denoise <- dfbm(count_mat = count_test_normalized,
                       quantiles = seq(0.1, 0.9, 0.1),
                       increment=j/10,
                       max_K=6,
                       ignore=0.05,
                       interpolate=FALSE,
                       ncores=4)
  count_test_denoised <- test_denoise$denoised_counts
  rownames(count_test_denoised) <- rownames(count_test_normalized)
  colnames(count_test_denoised) <- colnames(count_test_normalized)
  output[[sprintf("test_denoised_%d", j)]] <- count_test_denoised
  
  
}


saveRDS(output, sprintf("denoised_24/denoised_24_%d.rds", myseed))


