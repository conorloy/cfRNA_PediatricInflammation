suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

## Read in meta data and subset
mdf = read.csv("../0_DATA/metadata_all.csv")
    
## read in count matrix and subset
cnts = read.csv("../0_DATA/counts_all.csv",row.names=1)
cnts <- cnts[,mdf$Xsample_id]

## read in gene key
genekey = read.delim("../0_DATA/gencode.biotype.name.key.tsv")

## get all pairwise comparisons
all_groups = unique(mdf$inflam_cat)
# print(all_groups)
# cat("\n")
comparisons = apply(combn(all_groups,2),2,paste,collapse="<>")
# print(comparisons)
comparisons = comparisons[grepl("Healthy",comparisons)]
comparisons = comparisons[!grepl("ther",comparisons)]

SEEDS_ALL = c(1:10)

# comparisons = comparisons[1]
# SEEDS_ALL = c(1:2)

## run pipeline
COMP_VAR = "inflam_cat"
output_list <- list()

for (SEED in SEEDS_ALL){
    set.seed(SEED)
    for (COMP in comparisons){

        GROUPS = c(unlist(strsplit(COMP,"<>"))[1],
                    unlist(strsplit(COMP,"<>"))[2])

        samples <- mdf %>% 
            filter(.data[[COMP_VAR]] %in% GROUPS)
        rownames(samples) <- samples$Xsample_id

        ### subsample
        case_samps_sub <- samples %>% filter(.data[[COMP_VAR]] %in% GROUPS[1]) 

        nsub = 9
        if(nrow(case_samps_sub) <= nsub){nsub = nrow(case_samps_sub)}
        case_samps_sub <- case_samps_sub %>% sample_n(nsub)

        control_samps_sub <- samples %>% filter(.data[[COMP_VAR]] %in% GROUPS[2]) 
        samples <- rbind(case_samps_sub,control_samps_sub)
   
        counts <- cnts[,rownames(samples)]

        ##------------------------------------
        # Run DESeq
        dds <- DESeqDataSetFromMatrix(counts,
                                        colData = samples,
                                        design = formula(paste0("~",COMP_VAR,"+0")))
        dds <- DESeq(dds)

        ##------------------------------------
        # Get Results
        res <- results(dds,contrast = c(COMP_VAR,GROUPS[1], GROUPS[2])) 
        res_df <- res %>% data.frame()
        res_df <- merge(res_df,genekey, by.x=0, by.y="gene_id")

        output_list[[paste0(COMP,"_",SEED)]] = list("res"=res,"res_df"=res_df,"dds"=dds)
    }
}

## save RDS
output_list %>% saveRDS("./sub_output/pairwise_DAA.sub.9.rds")