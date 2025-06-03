
library(phyloseq)
library(ADAPT)
rm(list=ls())
data_12 <- readRDS("ECC/phyasv_visit12.rds")
data_24 <- readRDS("ECC/phyasv_visit24.rds")
# metadata_12 <- sample_data(data_12)

library(ADAPT)

da_result_12 <- adapt(data_12, cond.var="CaseEver")
da_result_24 <- adapt(data_24, cond.var="CaseEver")

result_table_12 <- da_result_12@details
result_table_24 <- da_result_24@details

write.csv(result_table_12, "ECC/DAA/DAA_result_yr1.csv",
          row.names=FALSE)

write.csv(result_table_24, "ECC/DAA/DAA_result_yr2.csv",
          row.names=FALSE)

marker_ASVs_12 <- result_table_12$Taxa[result_table_12$pval < 0.05]
marker_ASVs_24 <- result_table_24$Taxa[result_table_24$pval < 0.05]
marker_ASVs <- list(ASV_12=marker_ASVs_12, ASV_24=marker_ASVs_24)
saveRDS(marker_ASVs, "ECC/DAA/marker_ASVs.rds")

