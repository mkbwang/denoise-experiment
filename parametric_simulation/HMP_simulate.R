
library(simlite)
library(parallelly)
library(optparse)
folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment"

option_list <- list(make_option(c("-n", "--number"), type="integer", default=50, 
                                help="Number of samples for each phenotype [default=50]"),
                    make_option(c("-d", "--dispersion"), type="integer", default=5, 
                                help="Dispersion size parameter for negative binomial distribution [default=5]"),
                    make_option(c("-s", "--seed"), type="integer", default=1, 
                                help="seed [default=1]"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
nsample <- opt$number
dispersion <- opt$dispersion
myseed <- opt$seed

ncores <- availableCores()-1

zinb_parameters <- readRDS(file.path(folder, 
                                     "parametric_simulation/parameters/zinb_parameters.rds"))



# simulate from zero inflated negative binomial distribution
count_mat_list <- list()
mean_mat_list <- list()

set.seed(myseed)
for (j in 1:length(zinb_parameters)){
  
  selected_parameter <- zinb_parameters[[j]]
  nfeature <- nrow(selected_parameter)
  selected_parameter[, "size"] <- dispersion
  
  mean_values <- (1-selected_parameter[, "pi0"]) * selected_parameter[, "mu"]
  mean_values_mat <- replicate(nsample, mean_values)
  mean_mat_list[[j]] <- mean_values_mat
  
  copulas <- gausscopula(n=nsample, Sigma=diag(nrow=nfeature))
  count_matrix <- simcountmat(copula=copulas, params=selected_parameter, 
                              dist="zinb", ncores=ncores)
  count_mat_list[[j]] <- count_matrix

}

mean_mat <- round(do.call(cbind, mean_mat_list), digits=4)
count_mat <- do.call(cbind, count_mat_list)

rownames(count_mat) <- rownames(mean_mat) <- sprintf("Feature%d", seq(1, nfeature))
colnames(count_mat) <- colnames(mean_mat) <- sprintf("Sample%d", seq(1, nsample*6))

output_count_file <- sprintf("parametric_simulation/data/sim_count_nb_n%d_%d.csv", 
                             nsample, myseed)
output_mean_file <- sprintf("parametric_simulation/data/mean_nb_n%d.csv", 
                            nsample)

write.csv(count_mat, file.path(folder, output_count_file), 
          quote=FALSE)
write.csv(mean_mat, file.path(folder, output_mean_file), quote=FALSE)


