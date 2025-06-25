library(mbDenoise)
# source("ZIPNVA.R")
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

# nsample <- 50
# mymodel <- "zinb"
# dispersion <- 5
# myseed <- 10

input_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/data"
input_file <- file.path(input_folder, mymodel,
                        sprintf("sim_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))

output_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/mbDenoise"
output_file <- file.path(output_folder, mymodel,
                         sprintf("denoised_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))


input_data <- read.csv(input_file, row.names=1)
input_mat <- as.matrix(input_data) |> t()

output <- ZIPPCApn(X=input_mat, V=NULL, family="negative.binomial",
            trace=T, n.factors=2, maxit=100, parallel=T)



denoised_counts <- output$muz
rownames(denoised_counts) <- rownames(input_mat)
colnames(denoised_counts) <- colnames(input_mat)
denoised_counts_df <- as.data.frame(round(denoised_counts, 4))

write.csv(denoised_counts_df, file=output_file, quote=FALSE)


