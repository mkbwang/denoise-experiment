
rm(list=ls())
library(SAVER)
library(optparse)

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


input_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/data"
input_file <- file.path(input_folder, mymodel,
                        sprintf("sim_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))

output_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/SAVER"
output_file <- file.path(output_folder, mymodel,
                         sprintf("denoised_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))



library(data.table)
library(parallelly)
num_cores <- parallelly::availableCores()


input_data <- read.csv(input_file, row.names=1) |> as.matrix()

# input_matrix <- input_data[]
saver_smoothed_counts <- saver(x=input_data, size.factor = 1,
                        ncores = num_cores, estimates.only = TRUE)

saver_smoothed_counts <- data.frame(t(saver_smoothed_counts))

write.csv(saver_smoothed_counts, file=output_file, quote=FALSE)


