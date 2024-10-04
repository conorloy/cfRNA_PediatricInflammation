##------------------------------------
## DEPENDENCIES AND REFERENCE FILES
##------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(yaml))
suppressMessages(library(caret))
suppressMessages(library(glmnet))

##------------------------------------
## INPUTS
##------------------------------------

meta_train_path <- snakemake@input[['meta_train']]
meta_test_path <- snakemake@input[['meta_test']]
meta_val_path <- snakemake@input[['meta_val']]

counts_train_path <- snakemake@input[['counts_train_vst_sig']]
counts_test_path <- snakemake@input[['counts_test_vst_sig']]
counts_val_path <- snakemake@input[['counts_val_vst_sig']]


output_prefix <- snakemake@params[['output_prefix']]
model_functions <- snakemake@params[['model_functions']]
MODEL <- snakemake@params[['model']]
SAMP_GROUP <- snakemake@params[['SAMP_GROUP']]
group1_name = unlist(strsplit(SAMP_GROUP,"<>"))[1]
group2_name = unlist(strsplit(SAMP_GROUP,"<>"))[2]

id_column <- snakemake@params[['id_column']]
group_column <- snakemake@params[['group_column']]


## load functions and prepare data
source(model_functions)


## load metadata and counts
meta_data_train = read.csv(meta_train_path) %>% 
    mutate(group = .data[[group_column]]) %>% filter(group %in% c(group1_name,group2_name)) %>% 
    mutate(group = factor(group , levels = c(group1_name, group2_name)))

meta_data_test = read.csv(meta_test_path) %>% 
    mutate(group = .data[[group_column]]) %>% filter(group %in% c(group1_name,group2_name)) %>% 
    mutate(group = factor(group , levels = c(group1_name, group2_name)))

meta_data_val = read.csv(meta_val_path) %>% 
    mutate(group = .data[[group_column]]) %>% filter(group %in% c(group1_name,group2_name)) %>% 
    mutate(group = factor(group , levels = c(group1_name, group2_name)))


count_matrix_train_sig = read.delim(counts_train_path)
count_matrix_test_sig = read.delim(counts_test_path)
count_matrix_val_sig = read.delim(counts_val_path)

count_matrix_train_sig = count_matrix_train_sig[meta_data_train[[id_column]],]
count_matrix_test_sig = count_matrix_test_sig[meta_data_test[[id_column]],]
count_matrix_val_sig = count_matrix_val_sig[meta_data_val[[id_column]],]


# meta_data_train_all = read.csv(meta_train_path)
# meta_data_test_all = read.csv(meta_test_path)
# count_matrix_train_sig_all = read.delim(counts_train_path)
# count_matrix_test_sig_all = read.delim(counts_test_path)

# meta_data_all <- rbind(meta_data_train_all, meta_data_test_all)

# meta_data_train_all %>% nrow() %>% print()
# count_matrix_train_sig_all %>% nrow() %>% print()



##----------------------------------------------------------------------------
## ALL MODELS
##----------------------------------------------------------------------------
if(MODEL == "GLMNETLasso"){
    MODEL_FIT = GLMNETLasso_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "C5"){
    MODEL_FIT = C5_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "EXTRATREES"){
    MODEL_FIT = EXTRATREES_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "GLM"){
    MODEL_FIT = GLM_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}


if(MODEL == "GLMNETRidge"){
    MODEL_FIT = GLMNETRidge_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "KNN"){
    MODEL_FIT = KNN_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "LDA"){
    MODEL_FIT = LDA_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "NB"){
    MODEL_FIT = NB_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "NNET"){
    MODEL_FIT = NNET_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "PAM"){
    MODEL_FIT = PAM_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "RF"){
    MODEL_FIT = RF_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "RPART"){
    MODEL_FIT = RPART_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "SVMLIN"){
    MODEL_FIT = SVMLIN_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "SVMRAD"){
    MODEL_FIT = SVMRAD_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            meta_data_all, count_matrix_all_sig,
                            group1_name,group2_name,
                            "prob")
}


##----------------------------------------------------------------------------
## SAVE OUTPUTS
##----------------------------------------------------------------------------

MODEL_LIST = PERF_LIST[["MODEL_LIST"]]
train_outputs = PERF_LIST[["train_outputs"]]
test_outputs = PERF_LIST[["test_outputs"]]
val_outputs = PERF_LIST[["val_outputs"]]

MODEL_LIST %>% saveRDS(snakemake@output[['rds_output']])
train_outputs %>% write.table(snakemake@output[['train_stats']],sep="\t",row.names=T, quote=F)
test_outputs %>% write.table(snakemake@output[['test_stats']],sep="\t",row.names=T, quote=F)
val_outputs %>% write.table(snakemake@output[['val_stats']],sep="\t",row.names=T, quote=F)

quit()



# ##----------------------------------------------------------------------------
# ## GLMNET LASSO
# ##----------------------------------------------------------------------------
# output_subfolder <- "GLMNETLasso"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.glmnet = GLMNETLasso_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.glmnet

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT


# ##----------------------------------------------------------------------------
# ## GLMNET RIDGE
# ##----------------------------------------------------------------------------
# output_subfolder <- "GLMNETRidge"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.glmnet = GLMNETRidge_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.glmnet

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT



# ##----------------------------------------------------------------------------
# ## RF
# ##----------------------------------------------------------------------------

# output_subfolder <- "RF"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.rf = RF_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.rf

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## SVM Linear
# ##----------------------------------------------------------------------------
# output_subfolder <- "SVMLin"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.svmlin = SVMLIN_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.svmlin

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## NB
# ##----------------------------------------------------------------------------
# output_subfolder <- "NB"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.nb = suppressWarnings(NB_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.nb

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## LDA
# ##----------------------------------------------------------------------------
# output_subfolder <- "LDA"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.lda = suppressWarnings(LDA_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.lda
# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT



# ##----------------------------------------------------------------------------
# ## RPART NET
# ##----------------------------------------------------------------------------
# output_subfolder <- "RPART"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.rpart = RPART_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.rpart

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT


# ##----------------------------------------------------------------------------
# ## KNN
# ##----------------------------------------------------------------------------
# output_subfolder <- "KNN"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.knn = KNN_TRAIN(count_matrix_train_sig,meta_data_train)
# MODEL_FIT <- fit.knn

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## EXTRATREES_TRAIN
# ##----------------------------------------------------------------------------
# output_subfolder <- "EXTRATREES"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.extratrees = suppressWarnings(EXTRATREES_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.extratrees

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## NNET
# ##----------------------------------------------------------------------------
# output_subfolder <- "NNET"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.nnet = suppressWarnings(NNET_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.nnet

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT


# ##----------------------------------------------------------------------------
# ## SVM Radial
# ##----------------------------------------------------------------------------
# output_subfolder <- "SVMRAD"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.svmrad = suppressWarnings(SVMRAD_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.svmrad

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## PAM
# ##----------------------------------------------------------------------------
# output_subfolder <- "PAM"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.pam = suppressWarnings(PAM_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.pam

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT

# ##----------------------------------------------------------------------------
# ## C5
# ##----------------------------------------------------------------------------
# output_subfolder <- "C5"
# cat(paste0("    |--> ",output_subfolder,"\n"))
# fit.c5 = suppressWarnings(C5_TRAIN(count_matrix_train_sig,meta_data_train))
# MODEL_FIT <- fit.c5

# meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
#                            meta_data_train,meta_data_test,meta_data_val,
#                            group1_name,group2_name,
#                            "prob")

# all_models[[output_subfolder]] <- MODEL_FIT