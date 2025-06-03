

rm(list=ls())
library(DFBM)
library(BenchmarkDenoise)

marker_ASVs <- readRDS("ECC/DAA/marker_ASVs.rds")

data_24 <- readRDS("ECC/phyasv_visit24.rds")
metadata_24 <- sample_data(data_24)
count_24 <- otu_table(data_24) |> data.frame() |> as.matrix()

# select marker genes
count_24 <- count_24[, marker_ASVs$ASV_24]
target <- metadata_24$CaseEver
target_binary <- 1*(target == "Case")

# normalize
count_24_normalized <- normalize(count_mat=count_24, libsize=1e3)

write.csv(data.frame(count_24_normalized),
          "ECC/denoised/count_24_marker_original.csv",
          quote=FALSE)

# denoise
for (j in 4:9){
  print(j)
  dfbm_result <- dfbm(count_mat=count_24_normalized,
                        max_K=6,
                        increment=j/10,
                        ignore=0.1,
                        interpolate=FALSE,
                        ncores=4)
  
  thresholds <- dfbm_result$thresholds
  count_24_denoised <- dfbm_result$denoised_counts
  saveRDS(list(thresholds=thresholds, denoised_counts=count_24_denoised),
          sprintf("ECC/denoised/count_24_marker_denoised_%d.rds", j))

}


