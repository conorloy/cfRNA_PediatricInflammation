##---------------------------------------------------------------------------------------
## Measure performance of each model - train and test
##---------------------------------------------------------------------------------------

meaure_performance <- function(MODEL_FIT,
                               count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                               meta_data_train,meta_data_test,meta_data_val,
                                meta_data_all, count_matrix_all_sig,
                               group1_name,group2_name,
                               TYPE){

    suppressMessages(library(pROC))
    suppressMessages(library(vioplot))

    ##------------------------------------
    # TRAIN
    ##------------------------------------
        
    ### calculate classifier score
    tpred = predict(MODEL_FIT, newdata = count_matrix_train_sig, type=TYPE)
    meta_data_train$classifier_score <- as.numeric(tpred[,1])

    ### Get AUC
    roc_auc = roc(meta_data_train$group ~ meta_data_train$classifier_score, plot = FALSE, print.auc = FALSE)
    auc = roc_auc$auc

    list("tpred" = tpred, "meta_data_train" = meta_data_train, "roc_auc" = roc_auc, "auc" = auc) %>% saveRDS("./tmp.rds")

    ### Get Youden
    youden_threshold <- coords(roc_auc, x="best", input="threshold", best.method="youden")$threshold

    ### Make predicitons 
    meta_data_train <- meta_data_train %>%
        mutate(pred = ifelse(classifier_score > youden_threshold,group1_name,group2_name)) %>%
        mutate(pred = factor(pred, levels = c(group1_name, group2_name))) %>%
        mutate(group = factor(group, levels = c(group1_name, group2_name)))

    ### get stats
    results <- confusionMatrix(data = meta_data_train$pred, 
                                reference = meta_data_train$group)

    if(length(youden_threshold)>1){
        if(length(youden_threshold) == 2){youden_threshold = youden_threshold[1]}
        else{
            youden_threshold <- median(youden_threshold)
        }
    }

    t1 = data.frame(row.names= c("auc","youden"),
                    value = c(auc, youden_threshold))
    t2 <- rbind(as.matrix(results,what="classes"),
            as.matrix(results,what="overall")) %>% 
            data.frame() %>% 
            rename(.,value = `.`)

    train_outputs <- rbind(t1,t2)

    ##------------------------------------
    # TEST
    ##------------------------------------
    ### calculate classifier score
    tpred = predict(MODEL_FIT, newdata = count_matrix_test_sig, type=TYPE)
    meta_data_test$classifier_score <- as.numeric(tpred[,1])

    meta_data_test_eval = meta_data_test %>% filter(group %in% c(group1_name, group2_name))

    ### Get AUC
    roc_auc = roc(meta_data_test_eval$group ~ meta_data_test_eval$classifier_score, plot = FALSE, print.auc = FALSE)
    auc = roc_auc$auc

    ### Make predicitons 
    meta_data_test <- meta_data_test %>%
        mutate(pred = ifelse(classifier_score > youden_threshold,group1_name,group2_name)) %>%
        mutate(pred = factor(pred, levels = c(group1_name, group2_name))) %>%
        mutate(group = factor(group, levels = c(group1_name, group2_name)))

    ### get stats
    results <- confusionMatrix(data = meta_data_test$pred, 
                                reference = meta_data_test$group)

    t1 = data.frame(row.names= c("auc","youden"),
                    value = c(auc,"N/A"))
    t2 <- rbind(as.matrix(results,what="classes"),
            as.matrix(results,what="overall")) %>% 
            data.frame() %>% 
            rename(.,value = `.`)

    test_outputs <- rbind(t1,t2)

    ##------------------------------------
    # VALIDATE
    ##------------------------------------

    ### calculate classifier score
    tpred = predict(MODEL_FIT, newdata = count_matrix_val_sig, type=TYPE)
    meta_data_val$classifier_score <- as.numeric(tpred[,1])

    ### Get AUC
    roc_auc = roc(meta_data_val$group ~ meta_data_val$classifier_score, plot = FALSE, print.auc = FALSE)
    auc = roc_auc$auc

    ### Make predicitons 
    meta_data_val <- meta_data_val %>%
        mutate(pred = ifelse(classifier_score > youden_threshold,group1_name,group2_name)) %>%
        mutate(pred = factor(pred, levels = c(group1_name, group2_name))) %>%
        mutate(group = factor(group, levels = c(group1_name, group2_name)))

    ### get stats
    results <- confusionMatrix(data = meta_data_val$pred, 
                                reference = meta_data_val$group)

    t1 = data.frame(row.names= c("auc","youden"),
                    value = c(auc,"N/A"))
    t2 <- rbind(as.matrix(results,what="classes"),
            as.matrix(results,what="overall")) %>% 
            data.frame() %>% 
            rename(.,value = `.`)

    val_outputs <- rbind(t1,t2)


    ##------------------------------------
    # ALL
    ##------------------------------------

    # ### calculate classifier score on all samples
    # tpred = predict(MODEL_FIT, newdata = count_matrix_all_sig, type=TYPE)
    # meta_data_all$classifier_score <- as.numeric(tpred[,1])


    ##------------------------------------
    # SAVE
    ##------------------------------------

    ### all data as an RDS
    MODEL_LIST <- list("model" = MODEL_FIT,
                 "meta_data_train" = meta_data_train,
                 "meta_data_test" = meta_data_test,
                 "meta_data_val" = meta_data_val#,
                #  "meta_data_all" = meta_data_all
                 )

    
    return(list("MODEL_LIST" = MODEL_LIST, 
                "train_outputs" = train_outputs, 
                "test_outputs" = test_outputs, 
                "val_outputs" = val_outputs))
 
} 

##---------------------------------------------------------------------------------------
## Models
##---------------------------------------------------------------------------------------


##------------------------------------
## LDA
##------------------------------------

LDA_TRAIN <- function(counts,meta_data){
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    counts$y <- factor(meta_data$group)

    lda_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    preProcess = c("center","scale"),
                    method='lda', 
                    trControl=control)
    
    return(lda_fit)

}




##------------------------------------
## LOGISTIC REGRESSION
##------------------------------------

GLM_TRAIN <- function(counts,meta_data){
    
    set.seed(42)
    
    # counts <- round(counts)
    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    glm_fit <- train(y~., 
                    data=counts, 
                    method='glm',
                    preProcess = c("center","scale"),
                     weights = model_weights,
                    family = "quasibinomial")
    
    return(glm_fit)

}

##------------------------------------
## GLMNET
##------------------------------------

GLMNET_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    tunegrid <- expand.grid(alpha = c(0.25,0.5,0.75),
                           lambda = seq(0.00001, 1, length=10000))
    
    counts$y <- factor(meta_data$group)

    glmnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='glmnet',
                    preProcess = c("center","scale"),
                    trControl=control,
                    tuneGrid = tunegrid)
    
    return(glmnet_fit)
    
}

##------------------------------------
## GLMNETRidge
##------------------------------------

GLMNETRidge_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    tunegrid <- expand.grid(alpha = 0,
                           lambda = seq(0.00001, 1, length=10000))
    
    counts$y <- factor(meta_data$group)

    glmnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='glmnet',
                    preProcess = c("center","scale"),
                    trControl=control,
                    tuneGrid = tunegrid)
    
    return(glmnet_fit)
    
}

##------------------------------------
## GLMNETLasso
##------------------------------------

GLMNETLasso_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    tunegrid <- expand.grid(alpha = 1,
                           lambda = seq(0.00001, 1, length=10000))
    
    counts$y <- factor(meta_data$group)

    glmnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='glmnet',
                    preProcess = c("center","scale"),
                    trControl=control,
                    tuneGrid = tunegrid)
    
    return(glmnet_fit)
    
}

##------------------------------------ X
## RANDOM FOREST
##------------------------------------

RF_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)

    metric <- "Accuracy"
    mtry <- sqrt(ncol(counts)-1)
    tunegrid <- expand.grid(.mtry = (1:75))
    
    counts$y <- meta_data$group

    rf_fit <- train(y~., 
                          data=counts, 
                          method='rf', 
                          metric=metric, 
                        preProcess = c("center","scale"),
                          tuneGrid = tunegrid, 
                          trControl=control)
    
    return(rf_fit)

}


##------------------------------------ X
## SUPPORT VECTOR MACHINE LINEAR
##------------------------------------

SVMLIN_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5,classProbs = TRUE)
    
    tunegrid <- expand.grid(C = seq(0.00001, 2, length = 1000))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    svm_fit <- train(y~., 
                    data=counts, 
                    method='svmLinear', 
                    metric = metric,
                    trControl=control,
                     tuneGrid = tunegrid,
                    preProcess = c("center","scale"))
    
    return(svm_fit)

}


##------------------------------------ X
## Naive Bayes
##------------------------------------

NB_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(fL = 0:5,
                            usekernel = c(TRUE,FALSE),
                            adjust = seq(0, 5, by = 1)
                            )

    counts$y <- factor(meta_data$group)

    nb_fit <- train(y~., 
                    data=counts, 
                    method='nb', 
                    metric = metric,
                    preProcess = c("center","scale"),
                    trControl=control,
                   tuneGrid = tunegrid)
    
    return(nb_fit)

}


##------------------------------------
## KNN
##------------------------------------

KNN_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"
#     tunegrid <- expand.grid(fL = 0:5,
#                             usekernel = c(TRUE,FALSE),
#                             adjust = seq(0, 5, by = 1)
#                             )

    control <- trainControl(method="cv",number=5)
    
    counts$y <- factor(meta_data$group)

    knn_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='knn', 
                    trControl=control,
                    preProcess = c("center","scale"),
                    tuneLength = 100                      # test different K's
                    )
    
    return(knn_fit)

}



##------------------------------------
## RPART
##------------------------------------

RPART_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    rpart_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='rpart', 
#                     preProcess = c("center","scale"),
                    trControl=control,
                    weights = model_weights,
                    tuneLength = 20)
    
    return(rpart_fit)

}


##------------------------------------
## XGBTREE
##------------------------------------

XGBTREE_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    xgbtree_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='xgbTree', 
                    preProcess = c("center","scale"),
                    weights = model_weights,
                    trControl=control)
    
    return(xgbtree_fit)

}

##------------------------------------
## NNET
##------------------------------------

NNET_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                            size = c(3, 5, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85,
                                    90, 95, 100, 105, 110, 115, 120))

    
    counts$y <- factor(meta_data$group)
    
   model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    nnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='nnet', 
                    preProcess = c("center","scale"),
                    trControl=control,
                    tuneGrid = tunegrid,
                      weights = model_weights)
    
    return(nnet_fit)

}

##------------------------------------
## SVM RADIAL
##------------------------------------

SVMRAD_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5,classProbs = TRUE)
    
    tunegrid <- expand.grid(C = seq(0, 2, length = 1000),
                           sigma = c(.01, .015, 0.2))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    svmrad_fit <- train(y~., 
                    data=counts, 
                    method='svmRadial', 
                    metric = metric,
                    trControl=control,
                    preProcess = c("center","scale"),
                    tuneGrid = tunegrid)
    
    return(svmrad_fit)

}

##------------------------------------
## PAM
##------------------------------------

PAM_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    pam_fit <- train(y~., 
                    data=counts, 
                    method='pam', 
                    metric = metric,
                    trControl=control,
                    tuneLength = 20,
                    preProcess = c("center","scale"))
    
    return(pam_fit)

}

##------------------------------------
## ADABOOST
##------------------------------------

ADABOOST_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(nlter = seq(0, 2, length = 20),
                           method = c(.01, .015, 0.2))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    aba_fit <- train(y~., 
                    data = counts, 
                    method = 'adaboost', 
                    metric = metric,
                    trControl = control,
                    tuneLength = 30,
                    preProcess = c("center","scale"))
    
    return(aba_fit)

}

##------------------------------------
## C5.0
##------------------------------------

C5_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(winnow = c(TRUE,FALSE), 
                            trials=c(1,5,10,15,20),
                            model="tree" )
    
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    c5_fit <- train(y~., 
                    data = counts, 
                    method = 'C5.0', 
                    metric = metric,
                    trControl = control,
                    tuneGrid = tunegrid,
                    weigths = model_weights,
                    preProcess = c("center","scale"))
    
    return(c5_fit)

}


##------------------------------------ X
## Extra Trees
##------------------------------------

EXTRATREES_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5, 
                            classProbs = TRUE)
        
    n = ncol(counts)
    if(ncol(counts) > 5){
        
        tunegrid <- expand.grid(mtry = seq(5, ncol(counts), by = 5),
                           min.node.size = seq(2,10),
                            splitrule = "extratrees"
                           )}else{
        tunegrid <- expand.grid(mtry = seq(1, ncol(counts), by = 1),
                           min.node.size = seq(2,10),
                            splitrule = "extratrees"
                           )}
                        

    metric <- "Accuracy"
    
    counts$y <- meta_data$group

    rf_fit <- train(y~., 
                          data=counts, 
                          method='ranger', 
                          metric=metric, 
                          tuneGrid = tunegrid, 
                          trControl=control,
#                           preProc = c("center", "scale"),
                        na.action = na.pass
                          )
    
    return(rf_fit)
}

##------------------------------------ 
## GFS
##------------------------------------

GREEDY_FORWARD_SEARCH <- function(counts,meta_data){
    
    set.seed(42)

    
}



