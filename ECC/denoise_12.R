

rm(list=ls())
library(DFBM)
library(BenchmarkDenoise)

marker_ASVs <- readRDS("ECC/DAA/marker_ASVs.rds")

data_12 <- readRDS("ECC/phyasv_visit12.rds")
metadata_12 <- sample_data(data_12)
count_12 <- otu_table(data_12) |> data.frame() |> as.matrix()

# select marker genes
count_12 <- count_12[, marker_ASVs$ASV_12]
target <- metadata_12$CaseEver
target_binary <- 1*(target == "Case")

# normalize
count_12_normalized <- normalize(count_mat=count_12, libsize=1e3)

write.csv(data.frame(count_12_normalized),
          "ECC/denoised/count_12_marker_original.csv",
          quote=FALSE)


# denoise
for (j in 4:9){
  print(j)
  dfbm_result <- dfbm(count_mat=count_12_normalized,
                        max_K=6,
                        increment=j/10,
                        ignore=0.1,
                        interpolate=FALSE,
                        ncores=4)
  
  thresholds <- dfbm_result$thresholds
  count_12_denoised <- dfbm_result$denoised_counts
  saveRDS(list(thresholds=thresholds, denoised_counts=count_12_denoised),
          sprintf("ECC/denoised/count_12_marker_denoised_%d.rds", j))

}


