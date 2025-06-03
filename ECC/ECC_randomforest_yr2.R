

rm(list=ls())

marker_ASVs <- readRDS("experiment/ECC/marker_ASVs.rds")
data_24 <- readRDS("experiment/ECC/phyasv_visit24.rds")


metadata_24 <- sample_data(data_24)
count_24 <- otu_table(data_24) |> data.frame() |> as.matrix()

target <- metadata_24$CaseEver

count_24_normalized <- count_24
for (j in 1:nrow(count_24)){
    count_24_normalized[j, ] <- count_24[j, ] / sum(count_24[j, ]) * 1e4
}


prevalences <- colMeans(count_24_normalized > 0)
count_24_normalized <- count_24_normalized[, prevalences >= 0.05]



# DFBM denoising
library(DFBM)

# count_24_dfbm_full <- count_24_normalized
count_24_marker <- count_24_normalized[, marker_ASVs$ASV_24]


dfbm_result <- dfbm(count_mat=count_24_marker, max_K=6,
                                     ignore=0.1, interpolate=FALSE,
                                     ncores=4)

count_24_dfbm_marker <- dfbm_result$denoised_counts
colnames(count_24_dfbm_marker) <- colnames(count_24_marker)
rownames(count_24_dfbm_marker) <- rownames(count_24_marker)


# try random forest

library(caret)
library(randomForest)
library(pROC)



library(doParallel)
library(parallelly)


repeat_rf <- function(X, y, times=20){
    
    
    cl <- makeCluster(availableCores() - 1)  # use one fewer than available cores
    registerDoParallel(cl)

    train_roc <- rep(0, times)
    test_roc <- rep(0, times)


    control <- trainControl(method = "cv",
                            number = 5,
                            summaryFunction = twoClassSummary,
                            classProbs = TRUE,        # needed for ROC
                            verboseIter = FALSE)
    tuneGrid <- expand.grid(mtry = c(2,3,4,5,6))
    var_importance_mat <- matrix(0, nrow=times, ncol=ncol(X))
    for (k in 1:times){
        print(k)
        set.seed(k)
        trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
        trainData <- X[trainIndex, ]
        testData  <- X[-trainIndex, ]
        train_labels <- y[trainIndex]
        test_labels <- y[-trainIndex]

        model.rf <- train(x=trainData,
                          y=as.factor(train_labels),
                          method = "rf",
                          metric = "ROC",    # optimize for AUROC
                          tuneGrid = tuneGrid,
                          trControl = control)
        
        vimp_result <- varImp(model.rf)
        var_importance_mat[k, ] <- vimp_result$importance$Overall
        
        train_roc[k] <- max(model.rf$results$ROC)
        predProbs <- predict(model.rf, newdata = testData, type = "prob")
        roc_obj <- roc(response = test_labels, predictor = predProbs[, 1])
        auroc <- auc(roc_obj)
        if (auroc < 0.5) auroc <- 1-auroc
        test_roc[k] <- auroc

    }

    stopCluster(cl)

    return(list(Train=train_roc, Test=test_roc, Variable_importance=var_importance_mat))

}



accuracy_original <- repeat_rf(X=count_24_marker, y=target,
                               times=50)

accuracy_dfbm <- repeat_rf(X=count_24_dfbm_marker, y=target,
                               times=50)

wilcox.test(accuracy_original$Test, accuracy_dfbm$Test, paired=T)

mean(accuracy_original$Test < accuracy_dfbm$Test)
View(cbind(accuracy_original$Test, accuracy_dfbm$Test))

test_aucs <- data.frame(Original= accuracy_original$Test,
                        Denoised=accuracy_dfbm$Test)


ggplot(test_aucs, aes(x=Original, y=Denoised)) + geom_point()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue")+
  scale_y_continuous(limits=c(0.65, 0.95),
                     breaks=seq(0.65, 0.95, 0.1)) +
  scale_x_continuous(limits=c(0.65, 0.95),
                     breaks=seq(0.65, 0.95, 0.1))


# 
# test_aucs$Original <- accuracy_original$Test
# test_aucs$Denoised <- accuracy_dfbm$Test

# 
# test_aucs <- rbind(test_aucs_original,
#                    test_aucs_dfbm_full,
#                    test_aucs_dfbm_marker)

# test_aucs$Type <- factor(test_aucs$Type,
#                          levels=c("Original", "Full Denoise", "Marker Denoise"))
# 
# 
# ggplot(test_aucs, aes(x=Type, AUC)) + geom_boxplot() +
#     xlab("Input") + ylab("Test AUC") + scale_y_continuous(limits=c(0.5, 0.95),
#                                                           breaks=seq(0.5, 0.95, 0.1))

# 
# compare_original_full <- data.frame(Original=test_aucs_original$AUC,
#                                     FullDenoise = test_aucs_dfbm_full$AUC)
# 
# 
# compare_original_marker <- data.frame(Original=test_aucs_original$AUC,
#                                     MarkerDenoise = test_aucs_dfbm_marker$AUC)







ggplot(compare_original_marker, aes(x=Original, y=MarkerDenoise)) + geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue")+
    scale_y_continuous(limits=c(0.55, 0.95),
                       breaks=seq(0.5, 0.9, 0.1)) +
    scale_x_continuous(limits=c(0.55, 0.95),
                       breaks=seq(0.5, 0.95, 0.1))

# PCOA plot

library(vegan)
bray_dist <- vegdist(count_24_normalized, method="bray")
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_coords <- data.frame(pcoa_res$points)
colnames(pcoa_coords) <- c("PCOA1", "PCOA2")
pcoa_coords$CaseStatus <- target

ggplot(pcoa_coords, aes(x=PCOA1, y=PCOA2, color=CaseStatus)) +
    geom_point()

bray_dist_dfbm <- vegdist(count_24_dfbm_marker, method="bray")
pcoa_res_dfbm <- cmdscale(bray_dist_dfbm, k = 2, eig = TRUE)

pcoa_coords_dfbm <- data.frame(pcoa_res_dfbm$points)
colnames(pcoa_coords_dfbm) <- c("PCOA1", "PCOA2")
pcoa_coords_dfbm$CaseStatus <- target


ggplot(pcoa_coords_dfbm, aes(x=PCOA1, y=PCOA2, color=CaseStatus)) +
    geom_point()

