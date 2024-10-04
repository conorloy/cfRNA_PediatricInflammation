##------------------------------------
## DEPENDENCIES AND REFERENCE FILES
##------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(rlang))
suppressMessages(library(magrittr))
suppressMessages(library(pROC))
suppressMessages(library(caret))
suppressMessages(library(yaml))
set.seed(42)


meta_train_path <- snakemake@input[['meta_train']]
meta_test_path <- snakemake@input[['meta_test']]
meta_val_path <- snakemake@input[['meta_val']]

counts_train_path <- snakemake@input[['counts_train_vst_sig']]
counts_test_path <- snakemake@input[['counts_test_vst_sig']]
counts_val_path <- snakemake@input[['counts_val_vst_sig']]

deseq_res_path <- snakemake@input[['deseq_results']]

output_prefix <- snakemake@params[['output_prefix']]
group_column <- snakemake@params[['group_column']]
id_column <- snakemake@params[['id_column']]
model_functions <- snakemake@params[['model_functions']]
MODEL <- snakemake@params[['model']]
SAMP_GROUP <- snakemake@params[['SAMP_GROUP']]
group1_name <- SAMP_GROUP
group2_name <- "Other"


## load functions and prepare data
source(model_functions)

## load metadata and counts
meta_data_train = read.csv(meta_train_path) %>% 
    filter(inflam_cat != "Other")  %>%
    mutate(group = ifelse(.data[[group_column]] == group1_name, group1_name, group2_name )) %>% 
    mutate(group = factor(group, levels = c(group1_name, group2_name)))

meta_data_test = read.csv(meta_test_path) %>% 
    filter(inflam_cat != "Other")  %>%
    mutate(group = ifelse(.data[[group_column]] == group1_name, group1_name, group2_name )) %>% 
    mutate(group = factor(group, levels = c(group1_name, group2_name)))

meta_data_val = read.csv(meta_val_path) %>% 
    filter(inflam_cat != "Other")  %>%
    mutate(group = ifelse(.data[[group_column]] == group1_name, group1_name, group2_name )) %>% 
    mutate(group = factor(group, levels = c(group1_name, group2_name)))

count_matrix_train_sig = read.delim(counts_train_path)
count_matrix_test_sig = read.delim(counts_test_path)
count_matrix_val_sig = read.delim(counts_val_path)

count_matrix_train_sig = count_matrix_train_sig[meta_data_train[[id_column]],]
count_matrix_test_sig = count_matrix_test_sig[meta_data_test[[id_column]],]
count_matrix_val_sig = count_matrix_val_sig[meta_data_val[[id_column]],]

## load DESeq results
deseq_res = read.delim(deseq_res_path)

## gene reference list
gene.list_file = "/workdir/cfrna/references/human/hg38/gencode.biotype.name.key.tsv"
gene.list = read.delim(gene.list_file, col.names=c("gene_id","gene_name","type"))


#-------------------------
## Set Genes
# gene_dir <- deseq_res %>% data.frame() %>% 
#     rownames_to_column(.,var = "gene_id") %>% 
#     mutate(diff = ifelse(log2FoldChange > 0, "up","down")) %>%
#     select(gene_id, diff) %>% unique()

genes <- colnames(count_matrix_train_sig) 

##------------------------------------------------------------------------
## Function
##------------------------------------------------------------------------
three_layer_gfs <- function(top_genes, candidate_genes, df, cor_plot_dat, SELECTION){

  best_score = 0
  last_best_score = 0
  while_var = 0

  if(length(candidate_genes) ==0){return(top_genes)}

  #------------------------------------------------
  ## Execute GFS on Training Data  
  #------------------------------------------------
  suppressWarnings(suppressMessages({
  while(while_var != 1){

    gfs_output <- list()

    #------------------------------------------------
    # Loop through each candidate gene
    for (gene in candidate_genes) {

      # Create a temporary vector combining top_genes and the new gene
      temp_genes <- c(top_genes, gene)
      
      # Construct the formula and fit the model
      formula <- as.formula(paste("group ~", paste(temp_genes, collapse = " + "), "- 1"))
      model_train <- glm(formula, data = df, family = "binomial")
      
      # Make predictions
      probs <- predict(model_train, newdata = df, type = "response")
      
      # Calculate roc object
      roc_train <- roc(df$group, probs)
      coords_df <- data.frame(coords(roc_train, "all"))

      # Calculate Youden's index for each row
      coords_df$youden <- coords_df$sensitivity + coords_df$specificity - 1

      # Choose the row with the highest Youden's index
      optimal <- coords_df[which.max(coords_df$youden),]

      # Get the Stats at the best Threshold based on Youden Index
      auc_score =  auc(roc_train)
      spec_score =  optimal$specificity
      sens_score =  optimal$sensitivity

      # Caclulate the absolute value of the average correation
      #  -of all genes currently in the model with the last gene
      current_gene_data = cor_plot_dat[[gene]]
      top_genes_data = cor_plot_dat[, top_genes]
      cor_values = cor(current_gene_data, top_genes_data)
      average_cor = mean(cor_values)
      current_cor <- abs(average_cor)

      gfs_output[[gene]] <- c(current_cor,spec_score, sens_score, auc_score)      
      }
    #------------------------------------------------

    #------------------------------------------------
    # Find best gene - if any
    gfs_res <- do.call("rbind", gfs_output) %>% data.frame()
    colnames(gfs_res) <- c("correlation","spec_score","sens_score","auc_score")

    ## if Specificity is selection criteria
    if (SELECTION == "SPEC"){
      top_performer <- gfs_res %>%
        arrange(desc(spec_score),correlation) %>%
        head(1)

      best_score =top_performer$spec_score

    ## if AUC is selection criteria
    } else if (SELECTION == "AUC"){
      top_performer <- gfs_res %>% filter(spec_score > 0.9) %>%
      arrange(desc(auc_score),correlation) %>%
      head(1)

      best_score = top_performer$auc_score

    ## if Sensitivity is selection criteria
    } else if (SELECTION == "SENS"){
      top_performer <- gfs_res %>% filter(spec_score > 0.9) %>%
      arrange(desc(sens_score),correlation) %>%
      head(1)

      best_score <- top_performer$sens_score
    } else(print("PROBLEM WITH SELECTION CRITERIA"))

    # select the best gene
    best_gene = top_performer %>% rownames()

    if (nrow(top_performer) == 0){
      while_var <- 1
    } else if ((best_score - last_best_score) > 0.01) {
      # if yes, add to gene list and save specificity
      top_genes <- c(top_genes, best_gene)
      last_best_score <- best_score
    } else {
      # if not finish while loop
      while_var <- 1 
    }
    #------------------------------------------------
  
  whiles = 1
  }
  }))

  return(top_genes)
}

##------------------------------------------------------------------------
## Execute GFSu
##------------------------------------------------------------------------

#train data
norm_plot_dat_train = merge(count_matrix_train_sig, meta_data_train, by.x=0, by.y="sample_id") %>%
    rename(sample_id = `Row.names`)

# Prepare metadata (group as factor)
meta_data_train_cnts = norm_plot_dat_train             
meta_data_train_cnts$group <- ifelse(meta_data_train_cnts$group == group1_name, 1, 0)
meta_data_train_cnts$group <- as.factor(meta_data_train_cnts$group)

# Start with GBP5
starter_genes = c()

# data with counts
meta_and_counts = norm_plot_dat_train
    
# Run iterative GFS
best_spec <- three_layer_gfs(starter_genes, genes, meta_data_train_cnts, meta_and_counts, SELECTION="SPEC")
best_spec_auc <- three_layer_gfs(best_spec, genes[!(genes %in% best_spec)], meta_data_train_cnts, meta_and_counts, SELECTION="AUC")
best_spec_auc_sens <- three_layer_gfs(best_spec_auc, genes[!(genes %in% best_spec_auc)], meta_data_train_cnts, meta_and_counts, SELECTION="SENS")

print(best_spec_auc_sens)

# Construct the formula and fit the model
gene_names = as.character(best_spec_auc_sens)
formula <- as.formula(paste("group ~", paste(gene_names, collapse = "+"), "- 1"))
model_train <- glm(formula, data = meta_data_train_cnts, family = "binomial")

##------------------------------------------------------------------------
## Save outputs
##------------------------------------------------------------------------

output_subfolder <- "GFSu"
TYPE="response"
MODEL_FIT = model_train

suppressMessages(library(pROC))
suppressMessages(library(vioplot))

##------------------------------------
# TRAIN
##------------------------------------

### calculate classifier score
tpred = predict(MODEL_FIT, newdata = count_matrix_train_sig, type=TYPE)
meta_data_train$classifier_score <- as.numeric(tpred)

### Get AUC
roc_auc = roc(meta_data_train$group ~ meta_data_train$classifier_score, plot = FALSE, print.auc = FALSE)
auc = roc_auc$auc

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

t1 = data.frame(row.names= c("auc","youden"),
                value = c(auc, youden_threshold))
t2 <- rbind(as.matrix(results,what="classes"),
        as.matrix(results,what="overall")) %>% 
        data.frame() %>% 
        dplyr::rename(.,value = `.`)

train_outputs <- rbind(t1,t2)

##------------------------------------
# TEST
##------------------------------------

### calculate classifier score
tpred = predict(MODEL_FIT, newdata = count_matrix_test_sig, type=TYPE)
meta_data_test$classifier_score <- as.numeric(tpred)

### Get AUC
roc_auc = roc(meta_data_test$group ~ meta_data_test$classifier_score, plot = FALSE, print.auc = FALSE)
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
                value = c(auc,NA))
t2 <- rbind(as.matrix(results,what="classes"),
        as.matrix(results,what="overall")) %>% 
        data.frame() %>% 
        dplyr::rename(.,value = `.`)

test_outputs <- rbind(t1,t2)

##------------------------------------
# VALIDATE
##------------------------------------

### calculate classifier score
tpred = predict(MODEL_FIT, newdata = count_matrix_val_sig, type=TYPE)
meta_data_val$classifier_score <- as.numeric(tpred)

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
                value = c(auc,NA))
t2 <- rbind(as.matrix(results,what="classes"),
        as.matrix(results,what="overall")) %>% 
        data.frame() %>% 
        dplyr::rename(.,value = `.`)

val_outputs <- rbind(t1,t2)

##------------------------------------
# SAVE
##------------------------------------

### all data as an RDS
MODEL_LIST <- list("model" = MODEL_FIT,
                "meta_data_train" = meta_data_train,
                "meta_data_test" = meta_data_test,
                "meta_data_val" = meta_data_val,
                "genes" = best_spec_auc_sens)



MODEL_LIST %>% saveRDS(snakemake@output[['rds_output']])
train_outputs %>% write.table(snakemake@output[['train_stats']],sep="\t",row.names=T, quote=F)
test_outputs %>% write.table(snakemake@output[['test_stats']],sep="\t",row.names=T, quote=F)
val_outputs %>% write.table(snakemake@output[['val_stats']],sep="\t",row.names=T, quote=F)

quit()
