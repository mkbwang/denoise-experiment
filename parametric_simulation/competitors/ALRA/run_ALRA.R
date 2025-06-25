library(ALRA)
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

output_folder <- "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/ALRA"
output_file <- file.path(output_folder, mymodel,
                        sprintf("denoised_count_%s_n%d_d%d_%d.csv", mymodel, nsample, dispersion, myseed))


example_counts <- read.csv(input_file, row.names=1) |> t()

normalized_counts <- normalize_data(example_counts)

start.time <- Sys.time()
alra_result <- alra(normalized_counts)
end.time <- Sys.time()

alra_counts <- alra_result$A_norm_rank_k_cor_sc
alra_counts <- exp(alra_counts) - 1
alra_counts_df <- as.data.frame(round(alra_counts, 6))
rownames(alra_counts_df) <- rownames(example_counts)
colnames(alra_counts_df) <- colnames(example_counts)

write.csv(alra_counts_df, file=output_file, quote=FALSE)

