
folder <- "parametric_simulation/competitors/DFBM/zinb"


for(j in 1:100){
    filename <- sprintf("denoised_count_zinb_n50_d2_%d_cap.csv",
                        j)
    if (!file.exists(file.path(folder, filename))){
      print(filename)
    }
    
    filename <- sprintf("denoised_count_zinb_n50_d5_%d_cap.csv",
                        j)
    if (!file.exists(file.path(folder, filename))){
      print(filename)
    }
}

