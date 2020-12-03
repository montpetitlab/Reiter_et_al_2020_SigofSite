SEED = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 224, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 
        38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 
        56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 
        74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
        91, 92, 93, 94, 95, 96, 97, 98, 99, 100]

rule all:
    input: 
        expand('rna_seq/ava_seeds/ava_train70_optimal_rf_seed{seed}.RDS', seed = SEED),
        expand('rna_seq/site_seeds/site_train70_optimal_rf_seed{seed}.RDS', seed = SEED)


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
    script: "rna_seq/ava_seeds/scripts/random_forests_train70_ava.R"

#########################################
## RNA-seq site models
#########################################

rule variable_selection_site:
    input: 
        counts = "rna_seq/site/summarized_counts.RDS",
        pomona = "rna_seq/ava_seeds/pomona_install.txt"
    output:
        vita_rds = "rna_seq/site_seeds/site_vita_all_seed{seed}.RDS",
        vita_vars = "rna_seq/site_seeds/site_vita_all_vars_seed{seed}.txt",
        counts_filt = "rna_seq/site_seeds/site_vita_all_filt_seed{seed}.csv"
    conda: 'rna_seq/ava_seeds/env.yml'
    script: "rna_seq/site_seeds/scripts/variable_selection_site.R"

rule random_forests_train70_site:
    input:
        counts_filt = "rna_seq/site_seeds/site_vita_all_filt_seed{seed}.csv"
    output:
        rec_pars = "rna_seq/site_seeds/site_train70_rec_pars_seed{seed}.tsv",
        optimal_rf = 'rna_seq/site_seeds/site_train70_optimal_rf_seed{seed}.RDS',
        confusion_plt = 'rna_seq/site_seeds/site_train70_acc_seed{seed}.pdf'
    conda: "rna_seq/ava_seeds/env.yml"
    script: "rna_seq/site_seeds/scripts/random_forests_train70_site.R"
