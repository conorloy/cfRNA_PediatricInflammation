
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('output/seed_91_sep/1_feature_selection/seed_91_sep.counts.train.VST.tsv', 'output/seed_91_sep/1_feature_selection/seed_91_sep.counts.test.VST.tsv', 'output/seed_91_sep/1_feature_selection/seed_91_sep.counts.val.VST.tsv', 'input/seed_91_sep/seed_91_sep.metadata.train.tsv', 'input/seed_91_sep/seed_91_sep.metadata.test.tsv', 'input/seed_91_sep/seed_91_sep.metadata.val.tsv', 'output/seed_91_sep/1_feature_selection/seed_91_sep.deseqResults.tsv', "counts_train_vst" = 'output/seed_91_sep/1_feature_selection/seed_91_sep.counts.train.VST.tsv', "counts_test_vst" = 'output/seed_91_sep/1_feature_selection/seed_91_sep.counts.test.VST.tsv', "counts_val_vst" = 'output/seed_91_sep/1_feature_selection/seed_91_sep.counts.val.VST.tsv', "meta_train" = 'input/seed_91_sep/seed_91_sep.metadata.train.tsv', "meta_test" = 'input/seed_91_sep/seed_91_sep.metadata.test.tsv', "meta_val" = 'input/seed_91_sep/seed_91_sep.metadata.val.tsv', "deseq_results" = 'output/seed_91_sep/1_feature_selection/seed_91_sep.deseqResults.tsv'),
    output = list('output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.rds', 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.train.txt', 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.test.txt', 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.val.txt', "rds_output" = 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.rds', "train_stats" = 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.train.txt', "test_stats" = 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.test.txt', "val_stats" = 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES.val.txt'),
    params = list('output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES', 'tb', './workflow/scripts/func2_train-models.r', 'EXTRATREES', 'positive', 'negative', "output_prefix" = 'output/seed_91_sep/EXTRATREES/seed_91_sep.EXTRATREES', "group_column" = 'tb', "model_functions" = './workflow/scripts/func2_train-models.r', "model" = 'EXTRATREES', "group1_name" = 'positive', "group2_name" = 'negative'),
    wildcards = list('seed_91_sep', 'EXTRATREES', "sample" = 'seed_91_sep', "model" = 'EXTRATREES'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("MODELS" = c('C5', 'EXTRATREES', 'GLM', 'GLMNETLasso', 'GLMNETRidge', 'KNN', 'LDA', 'NB', 'NNET', 'PAM', 'RF', 'RPART', 'SVMLIN', 'SVMRAD'), "INPUT" = 'input/', "OUTPUT" = 'output/', "SAMPLES" = c('seed_1_sep', 'seed_2_sep', 'seed_3_sep', 'seed_4_sep', 'seed_5_sep', 'seed_6_sep', 'seed_7_sep', 'seed_8_sep', 'seed_9_sep', 'seed_10_sep', 'seed_11_sep', 'seed_12_sep', 'seed_13_sep', 'seed_14_sep', 'seed_15_sep', 'seed_16_sep', 'seed_17_sep', 'seed_18_sep', 'seed_19_sep', 'seed_20_sep', 'seed_21_sep', 'seed_22_sep', 'seed_23_sep', 'seed_24_sep', 'seed_25_sep', 'seed_26_sep', 'seed_27_sep', 'seed_28_sep', 'seed_29_sep', 'seed_30_sep', 'seed_31_sep', 'seed_32_sep', 'seed_33_sep', 'seed_34_sep', 'seed_35_sep', 'seed_36_sep', 'seed_37_sep', 'seed_38_sep', 'seed_39_sep', 'seed_40_sep', 'seed_41_sep', 'seed_42_sep', 'seed_43_sep', 'seed_44_sep', 'seed_45_sep', 'seed_46_sep', 'seed_47_sep', 'seed_48_sep', 'seed_49_sep', 'seed_50_sep', 'seed_51_sep', 'seed_52_sep', 'seed_53_sep', 'seed_54_sep', 'seed_55_sep', 'seed_56_sep', 'seed_57_sep', 'seed_58_sep', 'seed_59_sep', 'seed_60_sep', 'seed_61_sep', 'seed_62_sep', 'seed_63_sep', 'seed_64_sep', 'seed_65_sep', 'seed_66_sep', 'seed_67_sep', 'seed_68_sep', 'seed_69_sep', 'seed_70_sep', 'seed_71_sep', 'seed_72_sep', 'seed_73_sep', 'seed_74_sep', 'seed_75_sep', 'seed_76_sep', 'seed_77_sep', 'seed_78_sep', 'seed_79_sep', 'seed_80_sep', 'seed_81_sep', 'seed_82_sep', 'seed_83_sep', 'seed_84_sep', 'seed_85_sep', 'seed_86_sep', 'seed_87_sep', 'seed_88_sep', 'seed_89_sep', 'seed_90_sep', 'seed_91_sep', 'seed_92_sep', 'seed_93_sep', 'seed_94_sep', 'seed_95_sep', 'seed_96_sep', 'seed_97_sep', 'seed_98_sep', 'seed_99_sep', 'seed_100_sep')),
    rule = 'models',
    bench_iteration = as.numeric(NA),
    scriptdir = '/local/workdir/cjl332/tb-cfrna/SEED_SELECTION/1_all/workflow/rules/../scripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
##------------------------------------
## DEPENDENCIES AND REFERENCE FILES
##------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(yaml))
suppressMessages(library(caret))
suppressMessages(library(glmnet))

counts_train_path <- snakemake@input[['counts_train_vst']]
counts_test_path <- snakemake@input[['counts_test_vst']]
counts_val_path <- snakemake@input[['counts_val_vst']]
meta_train_path <- snakemake@input[['meta_train']]
meta_test_path <- snakemake@input[['meta_test']]
meta_val_path <- snakemake@input[['meta_val']]

deseq_res <- snakemake@input[['deseq_results']]

output_prefix <- snakemake@params[['output_prefix']]
group_column <- snakemake@params[['group_column']]
model_functions <- snakemake@params[['model_functions']]
MODEL <- snakemake@params[['model']]
group1_name <- snakemake@params[['group1_name']]
group2_name <- snakemake@params[['group2_name']]

## load functions and prepare data
source(model_functions)


## load metadata and counts
meta_data_train = read.delim(meta_train_path) %>% mutate(group = .data[[group_column]])
meta_data_test = read.delim(meta_test_path) %>% mutate(group = .data[[group_column]])
meta_data_val = read.delim(meta_val_path) %>% mutate(group = .data[[group_column]])

count_matrix_train_sig = read.delim(counts_train_path)
count_matrix_test_sig = read.delim(counts_test_path)
count_matrix_val_sig = read.delim(counts_val_path)

##----------------------------------------------------------------------------
## ALL MODELS
##----------------------------------------------------------------------------
if(MODEL == "C5"){
    MODEL_FIT = C5_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "EXTRATREES"){
    MODEL_FIT = EXTRATREES_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "GLM"){
    MODEL_FIT = GLM_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "GLMNETLasso"){
    MODEL_FIT = GLMNETLasso_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "GLMNETRidge"){
    MODEL_FIT = GLMNETRidge_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "KNN"){
    MODEL_FIT = KNN_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "LDA"){
    MODEL_FIT = LDA_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "NB"){
    MODEL_FIT = NB_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "NNET"){
    MODEL_FIT = NNET_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "PAM"){
    MODEL_FIT = PAM_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "RF"){
    MODEL_FIT = RF_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "RPART"){
    MODEL_FIT = RPART_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "SVMLIN"){
    MODEL_FIT = SVMLIN_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
                            group1_name,group2_name,
                            "prob")
}

if(MODEL == "SVMRAD"){
    MODEL_FIT = SVMRAD_TRAIN(count_matrix_train_sig,meta_data_train)
    PERF_LIST = meaure_performance(MODEL_FIT,
                            count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                            meta_data_train,meta_data_test,meta_data_val,
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

MODEL_LIST %>% saveRDS(paste0(output_prefix,".rds"))
train_outputs %>% write.table(paste0(output_prefix,".train.txt"),sep="\t",row.names=T, quote=F)
test_outputs %>% write.table(paste0(output_prefix,".test.txt"),sep="\t",row.names=T, quote=F)
val_outputs %>% write.table(paste0(output_prefix,".val.txt"),sep="\t",row.names=T, quote=F)

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