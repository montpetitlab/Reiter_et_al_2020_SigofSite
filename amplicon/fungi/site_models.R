library(dplyr)
library(readr)
library(tidyr)
library(ranger)
library(caret)
library(mlr)
library(tuneRanger)

set.seed(5)

fun_dada <- read_tsv("~/github/pn_amp/2019/outputs/ITS/dada2/ASV_counts.tsv") %>%
  select(-`16mockcomm1a`, -`16mockcomm2a`, -`16mockcomm3a`, -`16negctrl1a`,  
         -`16negctrl2a`, -`16negctrl3a`, -`17mockcomm1b`, -`17mockcomm2b`,
         -`17mockcomm3b`, -`17negctrl1b`, -`17negctrl2b`,  -`17negctrl3b`, 
         -`JFW_19-1-2`, -`JFW_19-19-2`, -`JFW_19-19-5`, -`JFW_19-2-2`, 
         -`JFW_19-20-5`, -`JFW_19-23-5`, -`JFW_19-24-2`, -`JFW_19-24-5`,
         -`negctrl1`, -`negctrl2`, -`negctrl3`, -`negctrl4`, -`19CH-50FNH`,
         -`19CH-50FOO`, -`19CH-50FPO`, -`19CH-70FOO`) %>%
  as.data.frame()

rownames(fun_dada) <- fun_dada$X1 # make ASVs rownames
fun_dada <- fun_dada[ , -1] # remove ASV column
fun_dada <- sweep(fun_dada, 2, colSums(fun_dada),`/`)
fun_dada <- t(fun_dada) # transpose
fun_dada <- as.data.frame(fun_dada)

tmp <- gsub("_[ABCD]{1}$", "", rownames(fun_dada))
tmp  <- gsub("201[679]{1}_", "", tmp)

table(tmp)
# missing: 2019 SNC2
# missing: 2019 RRV1
# missing: 2016 AV2
# missing: 2016 OR2

# SITE MODELS -------------------------------------------------------------
site_model <- function(counts, validation_year, rec_pars_file, optimal_rf_rds) {
  source("scripts/evaluate_model.R")
  # counts = fun_dada
  # validation_year = "2019"
  # rec_pars_file = "amplicon/fungi/site_train1617_rec_pars.tsv"
  # optimal_rf_rds = "amplicon/fungi/optimal_rf_site_train1617.RDS"
  ## separate years into train and validation sets
  counts_train <- !grepl(pattern = validation_year, x = rownames(counts))
  counts_train <- counts[counts_train, ]
  counts_val <- grepl(pattern = validation_year, x = rownames(counts))
  counts_val <- counts[counts_val, ]
  
  ## create vector with site information 
  counts_train$site <- gsub("_[ABCD]{1}$", "", rownames(counts_train))
  counts_train$site <- gsub("201[679]{1}_", "", counts_train$site)
  counts_val$site <- gsub("_[ABCD]{1}$", "", rownames(counts_val))
  counts_val$site <- gsub("201[679]{1}_", "", counts_val$site)
  
  ## tune/optimize model
  tmp <- counts_train
  colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
  task <- makeClassifTask(data = tmp, target = "site") # make an mlr task with counts17
  res <- tuneRanger(task, num.threads = 3)             # run tuning process
  write_tsv(res$recommended.pars, rec_pars_file)       # write model parameters to a file
  
  ## extract model parameters and use to build an optimal RF
  counts_train$site <- as.factor(counts_train$site)                # convert response to a factor
  optimal_rf <- ranger(
    dependent.variable.name = "site",
    mtry            = res$recommended.pars$mtry,
    num.trees        = 10000,
    data             = counts_train,
    sample.fraction = res$recommended.pars$sample.fraction,
    min.node.size   = res$recommended.pars$min.node.size,
    seed             = 5,
    importance       = 'permutation',
    local.importance = T
  )
  saveRDS(optimal_rf, file = optimal_rf_rds)
  # evaluate the accuracy of the model and generate a confusion matrix
  counts_train$site <- as.factor(counts_train$site)
  evaluate_model(optimal_ranger = optimal_rf,  data = counts_train, reference_class = counts_train$site, 
                 plt_title = "Training Performance") 
  
  # validation data
  counts_val$site <- factor(counts_val$site, levels = levels(counts_train$site))
  evaluate_model(optimal_ranger = optimal_rf, data = counts_val, reference_class = counts_val$site, 
                 plt_title = "Validation Performance") 
  
  return(optimal_rf)
}

# train on 16/17, validate on 19 
site_train1617 <- site_model(counts = fun_dada,
                             validation_year = "2019",
                             rec_pars_file = "amplicon/fungi/site_train1617_rec_pars.tsv",
                             optimal_rf_rds = "amplicon/fungi/optimal_rf_site_train1617.RDS")

# predicted
# observed AS1 AS2 AV1 AV2 CRN1 OR1 OR2 RRV1 RRV2 RRV3 SMV1 SMV2 SNC1 SNC2 SRH1
# AS1    4   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# AS2    0   4   0   0    0   0   0    0    0    0    0    0    0    0    0
# AV1    0   0   0   0    0   0   0    0    0    0    4    0    0    0    0
# AV2    0   0   0   0    0   0   0    0    0    0    0    3    0    0    0
# CRN1   0   0   0   0    1   0   0    0    0    0    0    0    0    3    0
# OR1    0   0   0   0    0   4   0    0    0    0    0    0    0    0    0
# OR2    0   0   0   0    0   3   1    0    0    0    0    0    0    0    0
# RRV1   0   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# RRV2   0   0   0   3    0   0   0    0    0    1    0    0    0    0    0
# RRV3   0   0   0   0    0   0   0    0    0    0    4    0    0    0    0
# SMV1   1   0   0   0    0   0   0    0    0    0    0    3    0    0    0
# SMV2   0   0   0   0    0   0   0    4    0    0    0    0    0    0    0
# SNC1   0   0   0   0    0   0   0    0    0    0    1    3    0    0    0
# SNC2   0   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# SRH1   0   0   0   0    0   0   0    0    0    0    1    0    0    0    3
# [1] "ACCURACY = 0.333333333333333"

# train on 16/19, validate on 17 
site_train1619 <- site_model(counts = fun_dada,
                             validation_year = "2017",
                             rec_pars_file = "amplicon/fungi/site_train1619_rec_pars.tsv",
                             optimal_rf_rds = "amplicon/fungi/optimal_rf_site_train1619.RDS")

# predicted
# observed AS1 AS2 AV1 AV2 CRN1 OR1 OR2 RRV1 RRV2 RRV3 SMV1 SMV2 SNC1 SNC2 SRH1
# AS1    1   0   0   0    3   0   0    0    0    0    0    0    0    0    0
# AS2    0   3   0   0    1   0   0    0    0    0    0    0    0    0    0
# AV1    0   0   0   0    4   0   0    0    0    0    0    0    0    0    0
# AV2    0   0   4   0    0   0   0    0    0    0    0    0    0    0    0
# CRN1   3   1   0   0    0   0   0    0    0    0    0    0    0    0    0
# OR1    0   0   0   0    0   4   0    0    0    0    0    0    0    0    0
# OR2    0   0   0   0    0   4   0    0    0    0    0    0    0    0    0
# RRV1   1   0   0   0    3   0   0    0    0    0    0    0    0    0    0
# RRV2   0   0   0   0    4   0   0    0    0    0    0    0    0    0    0
# RRV3   0   0   0   0    4   0   0    0    0    0    0    0    0    0    0
# SMV1   0   0   0   0    0   0   0    0    0    0    2    2    0    0    0
# SMV2   0   0   0   0    0   0   0    0    0    0    0    4    0    0    0
# SNC1   0   0   4   0    0   0   0    0    0    0    0    0    0    0    0
# SNC2   1   0   0   0    3   0   0    0    0    0    0    0    0    0    0
# SRH1   0   0   4   0    0   0   0    0    0    0    0    0    0    0    0
# [1] "ACCURACY = 0.233333333333333"

# train on 17/19, validate on 16 
site_train1719 <- site_model(counts = fun_dada,
                             validation_year = "2016",
                             rec_pars_file = "amplicon/fungi/site_train1719_rec_pars.tsv",
                             optimal_rf_rds = "amplicon/fungi/optimal_rf_site_train1719.RDS")


# predicted
# observed AS1 AS2 AV1 AV2 CRN1 OR1 OR2 RRV1 RRV2 RRV3 SMV1 SMV2 SNC1 SNC2 SRH1
# AS1    4   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# AS2    1   3   0   0    0   0   0    0    0    0    0    0    0    0    0
# AV1    0   0   4   0    0   0   0    0    0    0    0    0    0    0    0
# AV2    0   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# CRN1   0   0   0   0    2   0   0    0    1    0    1    0    0    0    0
# OR1    0   0   0   0    0   4   0    0    0    0    0    0    0    0    0
# OR2    0   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# RRV1   0   0   0   0    0   0   0    1    0    0    3    0    0    0    0
# RRV2   0   0   0   0    0   0   0    0    4    0    0    0    0    0    0
# RRV3   0   0   0   0    0   0   0    0    3    1    0    0    0    0    0
# SMV1   0   0   0   0    0   0   0    0    0    0    4    0    0    0    0
# SMV2   0   0   0   0    0   0   0    0    0    1    3    0    0    0    0
# SNC1   0   0   0   0    0   0   0    0    4    0    0    0    0    0    0
# SNC2   0   0   0   0    1   0   0    0    2    0    1    0    0    0    0
# SRH1   0   0   0   0    0   0   0    0    0    0    2    0    0    0    2
# [1] "ACCURACY = 0.557692307692308"

# train on 70% split, test on 30% -----------------------------------------
split_site_model <- function(counts, rec_pars_file, optimal_rf_rds) {
  source("scripts/evaluate_model.R")
  
  ## create vector with site information 
  site <- gsub("_[ABCD]{1}$", "", rownames(counts))
  site <- gsub("201[679]{1}_", "", site)
  
  ## separate years into train and validation sets
  ## semi-randomly create test set
  train_ind <- createDataPartition(site, p = .7, list = FALSE, times = 1)
  train <- counts[train_ind, ]
  test <- counts[-train_ind, ]
  train$site <- site[train_ind]
  test$site <- site[-train_ind]
  
  ## tune/optimize model
  tmp <- train
  colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
  task <- makeClassifTask(data = tmp, target = "site") # make an mlr task with counts17
  res <- tuneRanger(task, num.threads = 3)             # run tuning process
  write_tsv(res$recommended.pars, rec_pars_file)       # write model parameters to a file
  
  ## extract model parameters and use to build an optimal RF
  train$site <- as.factor(train$site)                # convert response to a factor
  optimal_rf <- ranger(
    dependent.variable.name = "site",
    mtry            = res$recommended.pars$mtry,
    num.trees       = 10000,
    data            = train,
    sample.fraction = res$recommended.pars$sample.fraction,
    min.node.size   = res$recommended.pars$min.node.size,
    seed             = 6,
    importance       = 'permutation',
    local.importance = T
  )
  saveRDS(optimal_rf, file = optimal_rf_rds)
  
  # evaluate the accuracy of the model and generate a confusion matrix
  evaluate_model(optimal_ranger = optimal_rf,  data = train, reference_class = train$site, 
                 plt_title = "Training Performance") 
  # validation data
  test$site <- as.factor(test$site)
  evaluate_model(optimal_ranger = optimal_rf, data = test, reference_class = test$site, 
                 plt_title = "Validation Performance") 
  
  return(optimal_rf)
}

site_train70 <- split_site_model(counts = fun_dada,
                                 rec_pars_file = "amplicon/fungi/site_train70_rec_pars.tsv",
                                 optimal_rf_rds = "amplicon/fungi/optimal_rf_site_train70.RDS")

# predicted
# observed AS1 AS2 AV1 AV2 CRN1 OR1 OR2 RRV1 RRV2 RRV3 SMV1 SMV2 SNC1 SNC2 SRH1
# AS1    3   0   0   0    0   0   0    0    0    0    0    0    0    0    0
# AS2    0   3   0   0    0   0   0    0    0    0    0    0    0    0    0
# AV1    0   0   3   0    0   0   0    0    0    0    0    0    0    0    0
# AV2    0   0   0   1    0   0   0    0    0    0    0    0    1    0    0
# CRN1   0   0   0   0    2   0   0    0    0    0    0    1    0    0    0
# OR1    0   0   0   0    0   3   0    0    0    0    0    0    0    0    0
# OR2    0   0   0   0    0   0   2    0    0    0    0    0    0    0    0
# RRV1   0   0   0   0    0   0   0    2    0    0    0    0    0    0    0
# RRV2   0   0   0   0    0   0   0    0    3    0    0    0    0    0    0
# RRV3   0   0   0   0    0   0   0    0    0    3    0    0    0    0    0
# SMV1   0   0   0   0    0   0   0    0    0    0    1    2    0    0    0
# SMV2   0   0   0   0    0   0   0    0    0    0    0    3    0    0    0
# SNC1   0   0   0   0    0   0   0    0    0    0    0    0    3    0    0
# SNC2   0   0   0   0    1   0   0    0    0    0    0    0    0    1    0
# SRH1   0   0   0   0    0   0   0    0    0    0    0    0    0    0    3
# [1] "ACCURACY = 0.878048780487805"