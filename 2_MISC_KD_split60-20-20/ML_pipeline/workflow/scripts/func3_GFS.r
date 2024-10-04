
#-------------------------
## Functions

## score calculator
calculate_counts <- function(counts_df, genes, dir){
    
    if (length(genes) == 0){
    
        scores <- data.frame(row.names = colnames(counts_df),rep(0,ncol(counts_df))) %>% set_colnames(paste0(dir,"_score"))
        colnames(scores)
        
        return(scores)
    }
    
    scores <- counts_df[rownames(counts_df) %in% genes,] %>% data.frame() %>% mutate_all(as.numeric) %>% 
                  mutate(log2( . + 1)) %>% summarise_all(mean) %>% t() %>% set_colnames(paste0(dir,"_score"))
    
    return(scores)
}



test_new_gene <- function(g,counts_df,combination,direction,metadata,sample_column,group_column){

    # add gene to gene_sub
    s = c(combination, g)
    s = s[!is.na(s)] # remvoe NAs

    # separate the genes into "up" and "down" groups
    up = direction[direction$diff == "up",] %>% filter(gene_id %in% s) %>% pull(gene_id)
    down = direction[direction$diff == "down",] %>% filter(gene_id %in% s) %>% pull(gene_id)

    # calculate scores for up and down genes separately
    up_counts = calculate_counts(counts_df,up,"up")
    down_counts = calculate_counts(counts_df,down,"down")

    # caculate final scores
    temp_df <- merge(up_counts, down_counts, by=0)
    temp_df <- merge(temp_df, metadata, by.x="Row.names", by.y=sample_column)
    temp_df$score <- temp_df$up_score - temp_df$down_score
    if (max(temp_df$down_score) ==0){temp_df$score = -temp_df$score}

    ## Calculate AUC using the score as a separtor
    auc = suppressMessages(auc(temp_df[[group_column]], as.numeric(temp_df$score)))
    return(cbind("combination" = paste(s, collapse="_"), "auc" = auc))
}




calc_metrics_gfs <- function(all_gfs,SUBSET){
    
    
    ## TRAIN
    train_gfs <- all_gfs %>% filter(set == SUBSET)

    auc = suppressMessages(auc(train_gfs$grp, as.numeric(train_gfs$score)))

    train_gfs_res <- confusionMatrix(factor(train_gfs$class_prediction), factor(train_gfs$grp))

    stats1 <- data.frame(value = as.matrix(train_gfs_res,what="overall"))%>% rownames_to_column("stat")
    stats2 <- data.frame(value=as.matrix(train_gfs_res,what="classes")) %>% rownames_to_column("stat")

    stats <- rbind(stats1,rbind(stats2,c("auc",auc))) %>% 
        column_to_rownames(.,"stat") %>%
        mutate(value = as.numeric(value))
    stats <- stats[c("auc","Accuracy","Sensitivity","Specificity"),,drop=F]
    stats$value <- round(stats$value,3)

    tab_metrics <- make_table(stats)

    return(tab_metrics)
    
}




calc_final_scores <- function(counts,metadata,best_up,best_down){
    up_cnts <- calculate_counts(counts, best_up, "up")
    down_cnts <- calculate_counts(counts, best_down, "down")

    # caculate final scores
    temp_df <- merge(up_cnts, down_cnts, by=0)
    temp_df <- merge(temp_df, metadata, by.x="Row.names", by.y=sample_column)
    temp_df$score <- temp_df$up_score - temp_df$down_score
    if (max(temp_df$up_score) ==0){temp_df$score = -temp_df$score}
    
    return(temp_df)
}

get_auc_youden <- function(score_df,group,youden=FALSE){
    
    score_df$group = score_df[group]
    
    if (youden == TRUE){
    
        png(file="tmp.png")
        rocobj1 <- plot.roc(score_df$group, score_df$score,
                            percent=TRUE,
                            ci=TRUE, boot.n=10000, boot.stratified=TRUE,
                            #	           show.thres=FALSE,
                            print.auc=TRUE,                  # compute AUC (of AUC by default)
                            #                   print.auc.pattern="%.2f",
                            print.thres="best",
                            print.thres.best.method="youden")
        youden_threshold <- coords(rocobj1, x="best", input="threshold", best.method="youden")$threshold
        dev.off()
        system("rm tmp.png")
        
    } else {youden_threshold = NA }
    
        roc_auc = roc(score_df$group ~ score_df$score, plot = FALSE, print.auc = TRUE)
    
    return(list("auc"=as.numeric(roc_auc$auc), "youden" = youden_threshold))
    
}


calc_metrics <- function(score_df,youden){
    score_df <- score_df %>% mutate(class_prediction = ifelse(score < youden,"other","MISC"))

    results <- confusionMatrix(factor(score_df$class_prediction), factor(score_df$group))

    all <- rbind(as.matrix(results,what="overall"),
                as.matrix(results,what="classes")) 

    accuracy = all[1,1] %>% as.numeric()
    sensitivity = all[8,1] %>% as.numeric()
    specificity = all[9,1] %>% as.numeric()
    
    return(list("accuracy" = accuracy,
                "sensitivity" = sensitivity,
                "specificity" = specificity))
}
      




calc_final_scores_fig <- function(counts,metadata,best_up,best_down){
    up_cnts <- calculate_counts(counts, best_up, "up")
    down_cnts <- calculate_counts(counts, best_down, "down")

    # caculate final scores
    temp_df <- merge(up_cnts, down_cnts, by=0)
    temp_df <- merge(temp_df, metadata, by.x="Row.names", by.y="sample_id")
    temp_df$score <- temp_df$up_score - temp_df$down_score
#     if (max(temp_df$up_score) ==0){temp_df$score = -temp_df$score}
    
    return(temp_df)
}

get_auc_youden_fig <- function(score_df,group,youden=FALSE){
    
    score_df$group = score_df[[group]]

    roc_auc = roc(score_df$group ~ score_df$score, plot = FALSE, print.auc = TRUE)

    if (youden == TRUE){
       youden_threshold <- coords(roc_auc, x="best", input="threshold", best.method="youden")$threshold
    } else {youden_threshold = NA }

    return(list("auc"=as.numeric(roc_auc$auc), "youden" = youden_threshold))
    
}