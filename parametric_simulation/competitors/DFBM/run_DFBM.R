
rm(list=ls())
library(DFBM)
library(optparse)
library(parallelly)
library(BenchmarkDenoise)

option_list <- list(make_option(c("-n", "--number"), type="integer", default=50, 
                                help="Number of samples for each phenotype [default=50]"),
                    make_option(c("-d", "--dispersion"), type="integer", default=5, 
                                help="Dispersion size parameter for negative binomial distribution [default=5]"),
                    make_option(c("-s", "--seed"), type="integer", default=1, 
                                help="seed [default=1]"),
                    make_option(c("-m", "--model"), type="character", default="zinb", 
                                help="seed [default=zinb]")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
nsample <- opt$number
dispersion <- opt$dispersion
mymodel <- opt$model
myseed <- opt$seed

# 
# nsample <- 50
# dispersion <- 1
# mymodel <- "zinb"
# myseed <- 50



input_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/data"
input_file <- file.path(input_folder, mymodel,
                        sprintf("sim_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))

output_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/DFBM"
output_file_nocap <- file.path(output_folder, mymodel,
                         sprintf("denoised_count_%s_n%d_d%d_%d_nocap.csv", mymodel, nsample, dispersion, myseed))
output_file_cap <- file.path(output_folder, mymodel,
                             sprintf("denoised_count_%s_n%d_d%d_%d_cap.csv", mymodel, nsample, dispersion, myseed))

input_data <- read.csv(input_file, row.names=1)
input_mat <- t(as.matrix(input_data))

# input_mat <- normalize(input_mat)
quantiles <- quantile(as.vector(input_mat), 0.9)



ncores <- parallelly::availableCores()

# denoise with no capped values
result_nocap <- dfbm(count_mat=input_mat,
               increment=0.8, max_K=6, ncores=ncores)
count_nocap <- result_nocap$denoised_counts

# denoise with capped values
result_cap <- dfbm(count_mat=input_mat,
                     increment=0.8, max_K=6, ncores=ncores,
                     cap=quantiles)
count_cap <- result_cap$denoised_counts


rownames(count_nocap) <- rownames(count_cap) <- rownames(input_mat) 
colnames(count_nocap) <- colnames(count_cap) <- colnames(input_mat)


write.csv(round(count_nocap, 5), output_file_nocap,
          quote=FALSE)

write.csv(round(count_cap, 5), output_file_cap,
          quote=FALSE)


