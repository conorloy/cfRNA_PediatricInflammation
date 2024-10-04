rule feature_selection_and_normalization:
    input:
        counts_train = config['counts_train_path'],
        counts_test = config['counts_test_path'],
        counts_val = config['counts_val_path'],
        counts_vst_train = config['counts_vst_train_path'],
        counts_vst_test = config['counts_vst_test_path'],
        counts_vst_val = config['counts_vst_val_path'],
        meta_train = config['meta_train_path'],
        meta_test = config['meta_test_path'],
        meta_val = config['meta_val_path']
    output:
        counts_train_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.train.VST.SIG.tsv',
        counts_test_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.test.VST.SIG.tsv',
        counts_val_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.val.VST.SIG.tsv',
        deseq_results =  OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.deseqResults.tsv'
    params:
        output_prefix= OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}',
        id_column = config['id_column'],
        group_column = config['group_column'],
        NUM_MARKERS= config['NUM_MARKERS'],
        OCCURANCE_RDS= config['OCCURANCE_RDS'],
        SAMP_GROUP= '{condition_pair}'
    script: 
        '../scripts/1_feature-selection.R'
    

rule models:
    input:
        counts_train_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.train.VST.SIG.tsv',
        counts_test_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.test.VST.SIG.tsv',
        counts_val_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.val.VST.SIG.tsv',
        meta_train = config['meta_train_path'],
        meta_test = config['meta_test_path'],
        meta_val = config['meta_val_path']
    output:
        rds_output = OUTPUT+'{condition_pair}/{model}/{condition_pair}.{model}.rds',
        train_stats = OUTPUT+'{condition_pair}/{model}/{condition_pair}.{model}.train.txt',
        test_stats = OUTPUT+'{condition_pair}/{model}/{condition_pair}.{model}.test.txt',
        val_stats = OUTPUT+'{condition_pair}/{model}/{condition_pair}.{model}.val.txt'
    params:
        output_prefix= OUTPUT+'{condition_pair}/{model}/{condition_pair}.{model}',
        id_column = config['id_column'],
        group_column = config['group_column'],
        model_functions= "./workflow/scripts/func2_train-models.r",
        model= "{model}",
        SAMP_GROUP= '{condition_pair}'
    script:
        '../scripts/2_train-models.R'


rule GFSu:
    input:
        counts_train_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.train.VST.SIG.tsv',
        counts_test_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.test.VST.SIG.tsv',
        counts_val_vst_sig = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.counts.val.VST.SIG.tsv',
        meta_train = config['meta_train_path'],
        meta_test = config['meta_test_path'],
        meta_val = config['meta_val_path'],
        deseq_results = OUTPUT + '{condition_pair}/1_feature_selection/{condition_pair}.deseqResults.tsv'
    output:
        rds_output = OUTPUT+'{condition_pair}/GFSu/GFSu_{condition_pair}.GFSu.rds',
        train_stats = OUTPUT+'{condition_pair}/GFSu/GFSu_{condition_pair}.GFSu.train.txt',
        test_stats = OUTPUT+'{condition_pair}/GFSu/GFSu_{condition_pair}.GFSu.test.txt',
        val_stats = OUTPUT+'{condition_pair}/GFSu/GFSu_{condition_pair}.GFSu.val.txt'
    params:
        output_prefix= OUTPUT+'{condition_pair}/GFSu/{condition_pair}.GFSu',
        id_column = config['id_column'],
        group_column = config['group_column'],
        model_functions= "./workflow/scripts/func2_train-models.r",
        SAMP_GROUP= '{condition_pair}'
    script:
        '../scripts/2_train-GFSu.R'
