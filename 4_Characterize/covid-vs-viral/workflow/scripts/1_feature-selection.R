##------------------------------------
## DEPENDENCIES AND REFERENCE FILES
##------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(yaml))
suppressMessages(library(DESeq2))
suppressMessages(library(pROC))

##------------------------------------
## INPUTS
##------------------------------------

meta_train_path <- snakemake@input[['meta_train']]
meta_test_path <- snakemake@input[['meta_test']]
meta_val_path <- snakemake@input[['meta_val']]

counts_train_path <- snakemake@input[['counts_train']]
counts_test_path <- snakemake@input[['counts_test']]
counts_val_path <- snakemake@input[['counts_val']]

counts_vst_train_path <- snakemake@input[['counts_vst_train']]
counts_vst_test_path <- snakemake@input[['counts_vst_test']]
counts_vst_val_path <- snakemake@input[['counts_vst_val']]


output_prefix <- snakemake@params[['output_prefix']]
id_column <- snakemake@params[['id_column']]
group_column <- snakemake@params[['group_column']]
NUM_OVERLAP <- snakemake@params[['NUM_OVERLAP']]
NUM_MARKERS <- snakemake@params[['NUM_MARKERS']]
OCCURANCE_RDS <- snakemake@params[['OCCURANCE_RDS']]

SAMP_GROUP <- snakemake@params[['SAMP_GROUP']]
g1 = unlist(strsplit(SAMP_GROUP,"<>"))[1]
g2 = unlist(strsplit(SAMP_GROUP,"<>"))[2]

get_auc <- function(GENE,CNTS,META,GRP,IDCOL){

    ## get counts for group of interest
    grp_cnts <- CNTS[,META[[IDCOL]] == GRP]

    ## get counts for all other groups
    other_cnts <- CNTS[,META[[IDCOL]] != GRP]

    ## get gene counts for group of interest
    grp_gene_cnts <- grp_cnts[GENE,]

    ## get gene counts for all other groups
    other_gene_cnts <- other_cnts[GENE,]

    ## get AUC
    auc <- suppressMessages(roc(cases = t(grp_gene_cnts),controls = t(other_gene_cnts),plot=FALSE)$auc)
    
    return(auc)
}

##------------------------------------
## LOAD DATA
##------------------------------------

## load metadata and counts
meta_data_train = read.csv(meta_train_path) %>% mutate(group = .data[[group_column]]) %>% filter(group %in% c(g1,g2))
meta_data_test = read.csv(meta_test_path) %>% mutate(group = .data[[group_column]]) %>% filter(group %in% c(g1,g2))
meta_data_val = read.csv(meta_val_path) %>% mutate(group = .data[[group_column]]) %>% filter(group %in% c(g1,g2))

## Counts raw
count_matrix_train = read.csv(counts_train_path,row.names=1)
count_matrix_test = read.csv(counts_test_path,row.names=1)
count_matrix_val = read.csv(counts_val_path,row.names=1)

count_matrix_train = count_matrix_train[,meta_data_train[[id_column]]]
count_matrix_test = count_matrix_test[,meta_data_test[[id_column]]]
count_matrix_val = count_matrix_val[,meta_data_val[[id_column]]]

## Counts VST
count_matrix_train_vst = read.csv(counts_vst_train_path,row.names=1)
count_matrix_test_vst = read.csv(counts_vst_test_path,row.names=1)
count_matrix_val_vst = read.csv(counts_vst_val_path,row.names=1)

count_matrix_train_vst = count_matrix_train_vst[,meta_data_train[[id_column]]]
count_matrix_test_vst = count_matrix_test_vst[,meta_data_test[[id_column]]]
count_matrix_val_vst = count_matrix_val_vst[,meta_data_val[[id_column]]]

## Get sig genes

##------------------------------------
## DESEQ
##------------------------------------

dds <- DESeqDataSetFromMatrix(count_matrix_train,
                                colData = meta_data_train,
                                design = formula(paste0("~group")))

pv_thresh = 0.05
fc_thresh = 1
bm_thresh = 150
au_thresh = 0

pv_thresh = 0.01
fc_thresh = 0.25
bm_thresh = 100
au_thresh = 0.65


dds <- DESeq(dds)
res <- results(dds) 
all_sig_df <- res %>% data.frame() %>% rownames_to_column(.,"gene_id")
sig_genes_df <-  all_sig_df %>%                    
                        filter(padj < pv_thresh) %>%                                             ## filter by p-values
                        filter(abs(log2FoldChange) > fc_thresh) %>%                                 ## filter by log2FC
                        filter(baseMean > bm_thresh)                                               ## filter by baseMean
                
## add auc
print(SAMP_GROUP)
sig_genes_df$auc <- apply(sig_genes_df,1,function(x) get_auc(x["gene_id"],count_matrix_train_vst,meta_data_train,g1,"group"))

## Filter and subset
sig_genes_df <- sig_genes_df %>% 
                        filter(auc >= au_thresh)  %>%                              ## TEMP !!!
                        arrange(desc(auc)) %>% head(NUM_MARKERS)                            ## filter by top N markers                        
sig_genes <- sig_genes_df %>% pull(gene_id)


## subset counts to sig genes
count_matrix_train_sig <- count_matrix_train_vst %>% 
        t() %>% data.frame() %>% select(all_of(sig_genes))
count_matrix_train_sig <- count_matrix_train_sig[meta_data_train[[id_column]],]

count_matrix_test_sig <- count_matrix_test_vst %>% 
        t() %>% data.frame() %>% select(all_of(sig_genes))  
count_matrix_test_sig <- count_matrix_test_sig[meta_data_test[[id_column]],]

count_matrix_val_sig <- count_matrix_val_vst %>% 
        t() %>% data.frame() %>% select(all_of(sig_genes))  
count_matrix_val_sig <- count_matrix_val_sig[meta_data_val[[id_column]],]

    
    
## Save files
count_matrix_train_sig %>% write.table(snakemake@output[['counts_train_vst_sig']], sep = "\t")
count_matrix_test_sig %>% write.table(snakemake@output[['counts_test_vst_sig']], sep = "\t")
count_matrix_val_sig %>% write.table(snakemake@output[['counts_val_vst_sig']], sep = "\t")
sig_genes_df %>% write.table(snakemake@output[['deseq_results']], sep = "\t")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# feature_filtering <- function(counts,
#                               CORR_THRESH = 0.75, 
#                               LOW_THRESH = 15
#                              ){
    
#     #####
#     ## input: 
#     # - counts = count matrix, normalized cpm (Columns are genes, rows are samples)
#     # - CORR_THRESH = max corerlation to not remove a gene
#     # - LOW_THRESH = minimum sum of all counts for that gene to be removed
#     ## output:
#     # - counts.cor.var.cnts = filtered count matrix
#     #####
    
#     suppressMessages(library(caret))
    
#     ## Convert to CPM for filtering
#     counts_cpm <- counts %>% edgeR::cpm() %>% t() %>% data.frame()
    
#     ## Select genes to remove based on counts
#     nzv <- nearZeroVar(counts_cpm)                                          # Variance

#     low.cnts <- c(1:ncol(counts_cpm))[colSums(counts_cpm) < LOW_THRESH]     # Low counts
    
#     thresh = dim(counts)[2]*0.75*0.5                                        # Threshold accross samples
# 	cut = counts %>% edgeR::cpm() %>% t() %>% colSums() %>% data.frame() %>% filter(. < thresh) %>% rownames() 
    
#     remove <- c(nzv,low.cnts,cut)

#     counts_sub <- counts[!(rownames(counts) %in% all_of(remove)),]
    
#     return(counts_sub)
    
# }

# ##------------------------------------
# ## PERFORM FEATURE SELECTION
# ##------------------------------------
# cat("--> PERFORMING FEATURE SELECTION\n")

# ##--------------------
# ## FUNCTIONS

# ## function to perform DESeq2
# run_deseq <- function(X_grp,
#                         mdf_grp,
#                         des){

#     # Contstruct DESeq Data Set
#     dds <- DESeqDataSetFromMatrix(X_grp,
#                                     colData = mdf_grp,
#                                     design = formula(des))
#     # DAA
#     dds <- DESeq(dds)
#     # Results
#     res <- results(dds,alpha=0.05,contrast = c("tb","positive", "negative")) %>% data.frame()

#     return(res)}

# ## function to get AUC
# get_auc <- function(gene,train_counts_cpm,train_meta){
#         gene_counts = train_counts_cpm[gene,,drop=FALSE] %>% t() 
#         meta_train_wGene <- merge(train_meta, gene_counts,by.x="sample_id",by.y=0)
#         roc_auc = roc(meta_train_wGene$tb ~ meta_train_wGene[[gene]], plot = FALSE, print.auc = FALSE, quiet=TRUE)
# return(as.numeric(roc_auc$auc))}

# ##--------------------
# ## DESeq2 on all samples

# count_matrix_train_filt <- feature_filtering(count_matrix_train)
# count_matrix_train_filt[is.na(count_matrix_train_filt)] <- 0
# res_all <- suppressMessages(run_deseq(as.matrix(count_matrix_train_filt), meta_data_train, "~tb + 0"))

# ##--------------------
# ## DESeq2 on separate cohorts

# ## If the cohort should be separated
# if(DAA_GRP == "sep"){

#     ## subset
#     mdf_trn_ghl  = meta_data_train %>% filter(cohort == "GHL")
#     mdf_trn_r2d2 = meta_data_train %>% filter(cohort == "R2D2")
#     mdf_trn_tbsq = meta_data_train %>% filter(cohort == "TBSQUARE")

#     X_trn_ghl  = count_matrix_train_filt[,mdf_trn_ghl$sample_id]
#     X_trn_r2d2 = count_matrix_train_filt[,mdf_trn_r2d2$sample_id]
#     X_trn_tbsq = count_matrix_train_filt[,mdf_trn_tbsq$sample_id]

#     ghl_res <- suppressMessages(run_deseq(X_trn_ghl, mdf_trn_ghl, "~tb + 0"))
#     r2d2_res <- suppressMessages(run_deseq(X_trn_r2d2, mdf_trn_r2d2, "~tb + 0"))
#     tbsq_res <- suppressMessages(run_deseq(X_trn_tbsq, mdf_trn_tbsq, "~tb + 0"))

#     ghl_sig <- ghl_res %>% filter(padj < 0.05) %>% rownames()
#     r2d2_sig <- r2d2_res %>% filter(padj < 0.05) %>% rownames()
#     tbsq_sig <- tbsq_res %>% filter(padj < 0.05) %>% rownames()

#     all_df <- data.frame(cohort = c( rep("ghl", length(ghl_sig)), rep("r2d2",length(r2d2_sig)), rep("tbsq",length(tbsq_sig)) ),
#                     genes = c(ghl_sig, r2d2_sig, tbsq_sig)
#                     )

#     sig_genes <- all_df %>% group_by(genes) %>% dplyr::count() %>% arrange(desc(n)) %>% filter(n>1) %>% pull(genes)

# }

# ##--------------------
# ## Get sig genes

# ## If the cohorts should be combined
# if(DAA_GRP == "all"){
#     sig_genes <- res_all %>% filter(padj < 0.05) %>% rownames()
# }

# ## if the length needs to be shortened
# if(length(sig_genes) > 100){

#     ## subset to sig genes
#     res_sig = res_all[sig_genes,] %>%
#         filter(baseMean > 50) %>% 
#         filter(padj < 0.05)

#     ### Calculate AUC per gene
#     res_sig$gene_auc <- as.numeric(lapply(rownames(res_sig),get_auc,
#                                                             count_matrix_train_vst,
#                                                             meta_data_train))
#     ## take top 100 by gene auc
#     sig_genes <- res_sig %>% 
#         arrange(desc(gene_auc)) %>%
#         head(100) %>%
#         rownames()

# }


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##------------------------------------
# ## SAVE
# ##------------------------------------

# ## subset to sig genes
# count_matrix_train_sig <- count_matrix_train_vst %>% 
#         t() %>% data.frame() %>% select(all_of(sig_genes))

# count_matrix_test_sig <- count_matrix_test_vst %>% 
#         t() %>% data.frame() %>% select(all_of(sig_genes))  

# count_matrix_val_sig <- count_matrix_val_vst %>% 
#         t() %>% data.frame() %>% select(all_of(sig_genes))  
    
    

# ## save RDS files

# tmp = data.frame("tmp" = "tmp")
# count_matrix_train_sig %>% write.table(paste0(output_prefix,".counts.train.VST.tsv"), sep = "\t")
# count_matrix_test_sig %>% write.table(paste0(output_prefix,".counts.test.VST.tsv"), sep = "\t")
# count_matrix_val_sig %>% write.table(paste0(output_prefix,".counts.val.VST.tsv"), sep = "\t")
# res_all %>% write.table(paste0(output_prefix,".deseqResults.tsv"), sep = "\t")

