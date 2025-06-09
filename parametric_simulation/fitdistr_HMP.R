
rm(list=ls())
library(BenchmarkDenoise)
library(simlite)


template <- readRDS("parametric_simulation/HMP_template.rds")

metadata <- template$metadata
count_mat <- template$count_mat



pois_list <- list()
zipois_list <- list()
nb_list <- list()
zinb_list <- list()
gamma_list <- list()
zigamma_list <- list()
lnorm_list <- list()
zilnorm_list <- list()


unique_types <- unique(metadata$HMP_BODY_SUBSITE)
nfeatures <- nrow(count_mat)

for (site_type in unique_types){
  
  print(site_type)
  
  count_mat_subset <- count_mat[, metadata$HMP_BODY_SUBSITE == site_type]
  prevalence <- rowMeans(count_mat_subset > 0)
  valid_features <- which(prevalence > 0.1) # do not fit models to features whose prevalence is smaller than 0.05
  
  ## poisson distribution
  print("Fit Poisson distribution")
  poisson_parameters <- matrix(0, nrow=nfeatures, ncol=1)
  poisson_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="pois",
                                 ncores=3)
  poisson_parameters[valid_features, ] <- poisson_fitted
  colnames(poisson_parameters) <- colnames(poisson_fitted)
  pois_list[[site_type]] <- poisson_parameters
  
  
  ## zero inflated poisson distribution
  print("Fit Zero Inflated Poisson distribution")
  zipois_parameters <- cbind(rep(1, nfeatures),
                             rep(0, nfeatures))
  zipois_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="zipois",
                                ncores=3)
  zipois_parameters[valid_features, ] <- zipois_fitted
  colnames(zipois_parameters) <- colnames(zipois_fitted)
  zipois_list[[site_type]] <- zipois_parameters
  
  
  ## negative binomial distribution
  print("Fit negative binomial distribution")
  nb_parameters <- cbind(matrix(0, nrow=nfeatures, ncol=2))
  nb_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="nb", ncores=3)
  nb_parameters[valid_features, ] <- nb_fitted
  colnames(nb_parameters) <- colnames(nb_fitted)
  nb_list[[site_type]] <- nb_parameters
  
  
  ## zero inflated nb distribution
  print("Fit zero-inflated negative binomial distribution")
  zinb_parameters <- cbind(rep(1, nfeatures),
                           matrix(0, nrow=nfeatures, ncol=2))
  zinb_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="zinb", ncores=3)
  zinb_parameters[valid_features, ] <- zinb_fitted
  colnames(zinb_parameters) <- colnames(zinb_fitted)
  zinb_list[[site_type]] <- zinb_parameters
  
  
  ## zero inflated gamma distribution
  print("Fit zero-inflated Gamma distribution")
  zigamma_parameters <- cbind(rep(1, nfeatures),
                              matrix(0, nrow=nfeatures, ncol=2))
  zigamma_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="zigamma", ncores=3)
  zigamma_parameters[valid_features, ] <- zigamma_fitted
  colnames(zigamma_parameters) <- colnames(zigamma_fitted)
  zigamma_list[[site_type]] <- zigamma_parameters
  
  
  ## zero inflated lognormal distribution
  print("Fit zero-inflated lognormal distribution")
  zilnorm_parameters <- cbind(rep(1, nfeatures),
                              matrix(0, nrow=nfeatures, ncol=2))
  zilnorm_fitted <- fitdistr_mat(count_mat_subset[valid_features, ], dist="zilnorm", ncores=3)
  zilnorm_parameters[valid_features, ] <- zilnorm_fitted
  colnames(zilnorm_parameters) <- colnames(zilnorm_fitted)
  zilnorm_list[[site_type]] <- zilnorm_parameters
  
}


saveRDS(pois_list, "parametric_simulation/parameters/pois_parameters.rds")
saveRDS(zipois_list, "parametric_simulation/parameters/zipois_parameters.rds")
saveRDS(nb_list, "parametric_simulation/parameters/nb_parameters.rds")
saveRDS(zinb_list, "parametric_simulation/parameters/zinb_parameters.rds")
saveRDS(zigamma_list, "parametric_simulation/parameters/zigamma_parameters.rds")
saveRDS(zilnorm_list, "parametric_simulation/parameters/zilnorm_parameters.rds")



