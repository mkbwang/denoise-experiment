

# first evaluate 24 month data

rm(list=ls())
library(BenchmarkDenoise)
library(phyloseq)
library(ggplot2)

phyasv_24 <- readRDS("ECC/phyasv_visit24.rds")

taxa_names <- tax_table(phyasv_24) |> as.data.frame()

marker_ASVs <- readRDS("ECC/DAA/marker_ASVs.rds")


metadata_24 <- sample_data(phyasv_24)

target <- metadata_24$CaseEver
target_binary <- 1*(target == "Case")

aucs_24_list <- list()
coefs_24_list <- list()
for (j in 1:100){
  
  result <- readRDS(sprintf("ECC/codalasso/performance_24_%d.rds", j))
  aucs_24_list[[j]] <- result$aucs_df
  coefs_24_list[[j]] <- result$coefs_mat
  
}

train_aucs <- lapply(aucs_24_list, function(x) x$train_auc)
train_aucs_matrix <- do.call(cbind, train_aucs)

test_aucs <- lapply(aucs_24_list, function(x) x$test_auc)
test_aucs_matrix <- do.call(cbind, test_aucs)


test_aucs_df <- t(test_aucs_matrix) |> as.data.frame()
colnames(test_aucs_df) <- rownames(aucs_24_list[[1]])

auc_plot <- ggplot(test_aucs_df, aes(x=Original, y=Denoise_7)) +
  geom_point(size=1, alpha=0.7) + 
  geom_abline(intercept=0, slope=1, linetype="dashed", color="blue")+
  scale_x_continuous(limits=c(0.5, 1),
                     breaks=seq(0.5,1, 0.1))+
  scale_y_continuous(limits=c(0.5, 1),
                     breaks=seq(0.5, 1, 0.1))+
  xlab("Original Test AUC") + ylab("Denoised Test AUC")



wilcox.test(test_aucs_df$Original, test_aucs_df$Denoise_7, paired=T)

# visualize MDS plot


count_24_normalized_original <- read.csv("ECC/denoised/count_24_marker_original.csv", 
                                         row.names=1) |> as.matrix()

dist_original <- dist_composition(count_mat=count_24_normalized_original,
                                  metric="bray", covariate=target)

mds_coordinates_original <- dist_original$mds_coordinates
mds_coordinates_original$Type <- "Before Denoising"


count_24_marker_denoised <- readRDS("ECC/denoised/count_24_marker_denoised_7.rds")
count_24_marker_denoised_mat <- count_24_marker_denoised$denoised_counts
rownames(count_24_marker_denoised_mat) <- rownames(count_24_normalized_original)
colnames(count_24_marker_denoised_mat) <- colnames(count_24_normalized_original)


num_entries <- rep(0, length(count_24_marker_denoised$thresholds))
for (j in 1:length(count_24_marker_denoised$thresholds)){
  num_entries[j] <- mean(count_24_normalized_original > count_24_marker_denoised$thresholds[j])
}


dist_denoised <- dist_composition(count_mat=count_24_marker_denoised_mat,
                                  metric="bray", covariate=target)


mds_coordinates_denoised <- dist_denoised$mds_coordinates
mds_coordinates_denoised$Type <- "After Denoising"

mds_coordinates_combined <- rbind(mds_coordinates_denoised, mds_coordinates_original)
colnames(mds_coordinates_combined)[3] <- "Status"

mds_coordinates_combined$Type <- factor(mds_coordinates_combined$Type,
                                        levels=c("Before Denoising", "After Denoising"))

mds_plot <- ggplot(mds_coordinates_combined, aes(x=V1, y=V2, color=Status)) + 
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab("PCOA1") + ylab("PCOA2") +
  facet_wrap(vars(Type))


coefs_original <- lapply(coefs_24_list, function(x) x["Original", ])
coefs_original <- do.call(rbind, coefs_original)
selection_probabilities_original <- colMeans(coefs_original != 0)

coefs_7 <- lapply(coefs_24_list, function(x) x["Denoise_7", ])
coefs_7 <- do.call(rbind, coefs_7)
selection_probabilities_7 <- colMeans(coefs_7 != 0) |> sort(decreasing=TRUE)


top_3_ASVs <- names(selection_probabilities_7[1:3])
ASV_names <- taxa_names[top_3_ASVs, ]

ASV_name_1 <- "Leptotrichia Wadei"
ASV_name_2 <- "Alloprevotella sp._HMT_914"
ASV_name_3 <- "Capnocytophaga sp._HMT_338"


count_24_topasv_original <- count_24_normalized_original[, top_3_ASVs] |> as.data.frame()
count_24_topasv_original$Status <- target

original_plot_1 <- ggplot(count_24_topasv_original, aes(x=ASV99, y=ASV101, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_1) + ylab(ASV_name_2)

original_plot_2 <- ggplot(count_24_topasv_original, aes(x=ASV99, y=ASV119, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_1) + ylab(ASV_name_3)

original_plot_3 <- ggplot(count_24_topasv_original, aes(x=ASV101, y=ASV119, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_2) + ylab(ASV_name_3)

library(patchwork)

original_plots <- (original_plot_1|original_plot_2|original_plot_3) + plot_layout(guides='collect')
print(original_plots)

count_24_topasv_denoised <- count_24_marker_denoised_mat[, top_3_ASVs] |> as.data.frame()
count_24_topasv_denoised$Status <- target

denoised_plot_1 <- ggplot(count_24_topasv_denoised, aes(x=ASV99, y=ASV101, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_1) + ylab(ASV_name_2)

denoised_plot_2 <- ggplot(count_24_topasv_denoised, aes(x=ASV99, y=ASV119, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_1) + ylab(ASV_name_3)

denoised_plot_3 <- ggplot(count_24_topasv_denoised, aes(x=ASV101, y=ASV119, color=Status))+
  geom_point(size=0.8, alpha=0.7) + 
  scale_color_manual(values=c("#A00000", "#0000A0"))+
  xlab(ASV_name_2) + ylab(ASV_name_3)  


denoised_plots <- (denoised_plot_1|denoised_plot_2|denoised_plot_3) + plot_layout(guides='collect')
print(denoised_plots)










