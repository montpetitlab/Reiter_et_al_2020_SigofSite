SEED = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

rule all:
    input: expand('rna_seq/ava_seeds/ava_train70_optimal_rf_seed{seed}.RDS', seed = SEED)


#########################################
## RNA-seq AVA models
#########################################

rule install_pomona:
    output:
        pomona = "rna_seq/ava_seeds/pomona_install.txt"
    conda: 'rna_seq/ava_seeds/env.yml'
    script: "rna_seq/ava_seeds/scripts/install_pomona.R"

rule variable_selection:
    input: 
        counts = "rna_seq/avas/summarized_counts_ava.RDS",
        pomona = "rna_seq/ava_seeds/pomona_install.txt"
    output:
        vita_rds = "rna_seq/ava_seeds/ava_vita_all_seed{seed}.RDS",
        vita_vars = "rna_seq/ava_seeds/ava_vita_all_vars_seed{seed}.txt",
        counts_filt = "rna_seq/ava_seeds/ava_vita_all_filt_seed{seed}.csv"
    conda: 'rna_seq/ava_seeds/env.yml'
    script: "rna_seq/ava_seeds/scripts/variable_selection.R"

rule random_forests_train70_ava:
    input:
        counts_filt = "rna_seq/ava_seeds/ava_vita_all_filt_seed{seed}.csv"
    output:
        rec_pars = "rna_seq/ava_seeds/ava_train70_rec_pars_seed{seed}.tsv",
        optimal_rf = 'rna_seq/ava_seeds/ava_train70_optimal_rf_seed{seed}.RDS',
        confusion_plt = 'rna_seq/ava_seeds/ava_train70_acc_seed{seed}.pdf'
    conda: "rna_seq/ava_seeds/env.yml"
    script: "rna_seq/ava_seeds/scripts/random_forets_train70_ava.R"


