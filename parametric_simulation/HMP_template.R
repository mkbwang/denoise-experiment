rm(list=ls())

library(scDesign2)
library(HMP16SData)
library(dplyr)

metadata <- V35() %>% table_one()
metadata_oral <- metadata %>% filter(`HMP Body Site` == "Oral")


V35_oral <- V35() %>% subset(select = HMP_BODY_SUBSITE %in%
                               c("Subgingival Plaque", "Supragingival Plaque", "Saliva",
                                 "Tongue Dorsum", "Throat", "Buccal Mucosa"))

oral_metadata <- colData(V35_oral) |> data.frame()
taxonomies <- rowData(V35_oral) |> data.frame()

counts_mat <- assay(V35_oral, "16SrRNA")
libsizes <- colSums(counts_mat)
prevalences <- rowMeans(counts_mat > 0)
taxa_filter <- prevalences > 0.2

# filter taxa based on prevalence
counts_mat<- counts_mat[taxa_filter, ]
libsizes <- colSums(counts_mat)

# filter sample based on library sizes
sample_filter <- libsizes > 5000 & libsizes < 20000
counts_mat <- counts_mat[, sample_filter]
oral_metadata <- oral_metadata[sample_filter, ]
prevalences <- rowMeans(counts_mat > 0)


which(counts_mat == max(counts_mat), arr.ind=TRUE)

# rank taxa by variance
count_var <- apply(counts_mat, 1, var)
variance_ranks <- rank(-count_var)
subset_otus <- which(variance_ranks <= 500)



# only select the most variable OTUs
counts_mat_subset <- counts_mat[subset_otus, ]

# reorder sample rows and columns
oral_metadata <- oral_metadata %>% arrange(HMP_BODY_SUBSITE)
counts_mat_subset <- counts_mat_subset[, rownames(oral_metadata)]


mapping <- c("Buccal Mucosa"="Buccal_mucosa", "Saliva"="Saliva", 
             "Subgingival Plaque"="Subgingival_Plaque", "Supragingival Plaque"="Supragingival_Plaque",
             "Throat"="Throat", "Tongue Dorsum" ="Tongue_Dorsum")
oral_metadata$HMP_BODY_SUBSITE <- mapping[oral_metadata$HMP_BODY_SUBSITE]


# simulate from this template using scDesign
colnames(counts_mat_subset) <- oral_metadata$HMP_BODY_SUBSITE
cell_type_prop <- table(oral_metadata$HMP_BODY_SUBSITE) / nrow(oral_metadata)

# visualize the template real data count matrix
log_real_counts <- log10(1+counts_mat_subset)

library(BenchmarkDenoise)
library(ComplexHeatmap)
library(circlize)

colors <- c("#1b4f72", "#4a7c59", "#7e1e9c", "#663399", "#d35400", "#444444")
names(colors) <- unique(oral_metadata$HMP_BODY_SUBSITE)

col_annotation <- HeatmapAnnotation(
  df = data.frame(Category=oral_metadata$HMP_BODY_SUBSITE),
  col=list(Category=colors)
)
col_fun <- colorRamp2(c(0, 2), c("#FFFFFF", "#AA0000"))
save_heatmap(X=log_real_counts,
             entry_name="Log Count",
             filename="scDesign/real_data_template.png",
             rowannot=NULL,
             colannot=col_annotation,
             cmap=col_fun,
             width=18, height=12,
             legend=T)


template <- list(metadata=oral_metadata,
                 count_mat=counts_mat_subset)

saveRDS(template, "parametric_simulation/HMP_template.rds")

# the rest are archived codes for simulation based on scDesign2
# unique_oral_sites <- unique(oral_metadata$HMP_BODY_SUBSITE)
# poisson_fit <- fit_model_scDesign2(data_mat = counts_mat_subset,
#                                    cell_type_sel=unique_oral_sites,
#                                    marginal="poisson",
#                                    zp_cutoff=1,
#                                    ncores=4)
# 
# 
# 
# poisson_sim <- simulate_count_scDesign2(model_params=poisson_fit,
#                                         n_cell_new = 200,
#                                         sim_method="copula",
#                                         cell_type_prop = cell_type_prop)
# 
# saveRDS(poisson_sim, file="experiment/scDesign/HMP/poisson_sim.rds")
# 
# 
# NB_fit <- fit_model_scDesign2(data_mat = counts_mat,
#                               cell_type_sel=c("Plaque", "Tongue", "Throat"),
#                               marginal="nb",
#                               zp_cutoff=1,
#                               ncores=2)
# 
# NB_sim <- simulate_count_scDesign2(model_params=NB_fit,
#                                         n_cell_new = 200,
#                                         sim_method="copula",
#                                         cell_type_prop = cell_type_prop)
# 
# saveRDS(NB_sim, file="experiment/scDesign/HMP/NB_sim.rds")
# 
# 
# ZINB_fit <- fit_model_scDesign2(data_mat = counts_mat,
#                                 cell_type_sel=c("Plaque", "Tongue", "Throat"),
#                                 marginal="zinb",
#                                 zp_cutoff=1,
#                                 ncores=2)
# 
# 
# ZINB_sim <- simulate_count_scDesign2(model_params=ZINB_fit,
#                                    n_cell_new = 200,
#                                    sim_method="copula",
#                                    cell_type_prop = cell_type_prop)
# 
# saveRDS(ZINB_sim, file="experiment/scDesign/HMP/ZINB_sim.rds")
# 
# 
# full_fit <- fit_model_scDesign2(data_mat = counts_mat,
#                                 cell_type_sel=c("Plaque", "Tongue", "Throat"),
#                                 marginal="auto_choose",
#                                 zp_cutoff=1,
#                                 ncores=2)
# 
# full_sim <- simulate_count_scDesign2(model_params=ZINB_fit,
#                                      n_cell_new = 200,
#                                      sim_method="copula",
#                                      cell_type_prop = cell_type_prop)
# 
# saveRDS(full_sim, file="experiment/scDesign/HMP/full_sim.rds")
# 
# 
# saveRDS(list(poisson=poisson_fit, NB=NB_fit,
#              ZINB=ZINB_fit, full=full_fit),
#         file="experiment/scDesign/HMP/fitted_parameters.rds")



