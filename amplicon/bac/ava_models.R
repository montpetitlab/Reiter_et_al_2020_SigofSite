library(dplyr)
library(readr)
library(tidyr)
library(ranger)
library(caret)
library(mlr)
library(tuneRanger)
library(stringr)
library(tibble)
set.seed(1)

bac_tax_contam <- read_tsv("~/github/pn_amp/2019/outputs/16s/dada2/ASV_taxonomy.tsv") %>%
  filter(Family == "Mitochondria" | Order == "Chloroplast")

bac_dada <- read_tsv("~/github/pn_amp/2019/outputs/16s/dada2/ASV_counts.tsv") %>%
  select(-`16mockcomm1a`, -`16mockcomm2a`, -`16mockcomm3a`, -`16negctrl1a`,  
         -`16negctrl2a`, -`16negctrl3a`, -`17mockcomm1b`, -`17mockcomm2b`,
         -`17mockcomm3b`, -`17negctrl1b`, -`17negctrl2b`,  -`17negctrl3b`, 
         -`JFW_19-1-2`, -`JFW_19-19-2`, -`JFW_19-19-5`, -`JFW_19-2-2`, 
         -`JFW_19-20-5`, -`JFW_19-23-5`, -`JFW_19-24-2`, -`JFW_19-24-5`,
         -`negctrl2`, -`negctrl3`, -`negctrl4`, -`19CH-50FNH`,
         -`19CH-50FOO`, -`19CH-50FPO`, -`19CH-70FOO`) %>%
  filter(!X1 %in% bac_tax_contam$X1) %>% # remove cholorplast and mitochondria
  as.data.frame() %>%
  column_to_rownames("X1")

bac_dada <- sweep(bac_dada, 2, colSums(bac_dada),`/`)
bac_dada <- t(bac_dada) # transpose
rownames(bac_dada) <- gsub("PN", "", rownames(bac_dada))
rownames(bac_dada) <- gsub("PMS", "MSA", rownames(bac_dada))
bac_dada <- as.data.frame(bac_dada)
rownames(bac_dada)
tmp <-gsub("[1-3]_[ABCD]{1}$", "", rownames(bac_dada))
tmp  <- gsub("201[679]_", "", tmp)

table(tmp)

# ava MODELS -------------------------------------------------------------
ava_model <- function(counts, validation_year, rec_pars_file, optimal_rf_rds) {
  source("scripts/evaluate_model.R")
  ## separate years into train and validation sets
  counts_train <- !grepl(pattern = validation_year, x = rownames(counts))
  counts_train <- counts[counts_train, ]
  counts_val <- grepl(pattern = validation_year, x = rownames(counts))
  counts_val <- counts[counts_val, ]
  
  ## create vector with ava information 
  counts_train$ava <- gsub("[1-3]_[ABCD]{1}$", "", rownames(counts_train))
  counts_train$ava <- gsub("201[679]_", "", counts_train$ava)
  counts_val$ava <- gsub("[1-3]_[ABCD]{1}$", "", rownames(counts_val))
  counts_val$ava <- gsub("201[679]_", "", counts_val$ava)
  
  
  ## tune/optimize model
  tmp <- counts_train
  colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
  task <- makeClassifTask(data = tmp, target = "ava")  # make an mlr task with counts17
  res <- tuneRanger(task, num.threads = 3)             # run tuning process
  write_tsv(res$recommended.pars, rec_pars_file)       # write model parameters to a file
  
  ## extract model parameters and use to build an optimal RF
  counts_train$ava <- as.factor(counts_train$ava)      # convert response to a factor
  optimal_rf <- ranger(
    dependent.variable.name = "ava",
    mtry            = res$recommended.pars$mtry,
    num.trees        = 10000,
    data             = counts_train,
    sample.fraction = res$recommended.pars$sample.fraction,
    min.node.size   = res$recommended.pars$min.node.size,
    seed             = 1,
    importance       = 'permutation',
    local.importance = T
    )
  saveRDS(optimal_rf, file = optimal_rf_rds)
  
  # evaluate the accuracy of the model and generate a confusion matrix
  evaluate_model(optimal_ranger = optimal_rf,  data = counts_train, reference_class = counts_train$ava, 
                 plt_title = "Training Performance") 
  # validation data
  counts_val$ava <- as.factor(counts_val$ava)
  evaluate_model(optimal_ranger = optimal_rf, data = counts_val, reference_class = counts_val$ava, 
                 plt_title = "Validation Performance") 
  
  return(optimal_rf)
}

# train on 16/17, validate on 19 
ava_train1617 <- ava_model(counts = bac_dada,
                           validation_year = "19",
                           rec_pars_file = "amplicon/bac/ava_train1617_rec_pars.tsv",
                           optimal_rf_rds = "amplicon/bac/optimal_rf_ava_train1617.RDS")

# predicted
# observed AS AV CRN OR RRV SMV SNC SRH
# AS   6  0   0  0   0   2   0   0
# AV   0  3   0  0   1   1   0   2
# CRN  0  0   0  0   0   3   1   0
# OR   0  0   0  8   0   0   0   0
# RRV  0  0   0  0   7   1   0   0
# SMV  0  0   0  0   0   8   0   0
# SNC  0  0   0  0   0   2   2   0
# SRH  0  0   0  0   0   0   0   4
# [1] "ACCURACY = 0.745098039215686"

# train on 16/19, validate on 17 
ava_train1619 <- ava_model(counts = bac_dada,
                           validation_year = "17",
                           rec_pars_file = "amplicon/bac/ava_train1619_rec_pars.tsv",
                           optimal_rf_rds = "amplicon/bac/optimal_rf_ava_train1619.RDS")

# predicted
# observed AS AV CRN OR RRV SMV SNC SRH
# AS   7  1   0  0   0   0   0   0
# AV   0  4   0  0   4   0   0   0
# CRN  0  0   4  0   0   0   0   0
# OR   0  0   0  8   0   0   0   0
# RRV  0  0   0  0  12   0   0   0
# SMV  0  4   0  0   0   3   0   1
# SNC  0  0   1  0   4   0   3   0
# SRH  0  4   0  0   0   0   0   0
# [1] "ACCURACY = 0.683333333333333"

# train on 17/19, validate on 16 
ava_train1719 <- ava_model(counts = bac_dada,
                           validation_year = "16",
                           rec_pars_file = "amplicon/bac/ava_train1719_rec_pars.tsv",
                           optimal_rf_rds = "amplicon/bac/optimal_rf_ava_train1719.RDS")
# predicted
# observed AS AV CRN OR RRV SMV SNC SRH
# AS   8  0   0  0   0   0   0   0
# AV   0  4   0  0   0   0   0   0
# CRN  0  0   0  0   0   0   4   0
# OR   0  0   0  3   1   0   0   0
# RRV  0  1   0  0  11   0   0   0
# SMV  3  0   0  0   0   4   0   1
# SNC  0  0   0  0   4   0   4   0
# SRH  2  0   0  0   0   0   0   2
# [1] "ACCURACY = 0.692307692307692"


# train on 70% split, test on 30% -----------------------------------------
split_ava_model <- function(counts, rec_pars_file, optimal_rf_rds) {
  source("scripts/evaluate_model.R")
  
  ## create vector with ava information 
  ava <- gsub("[1-3]_[ABCD]{1}$", "", rownames(counts))
  ava <- gsub("201[679]_", "", ava)
  
  ava <- ifelse(str_detect(string = ava, pattern = "P5A"), "AS", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "MSA"), "AS", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "ANA"), "SNC", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "CLD"), "SNC", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "CHW"), "CRN", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "NEL"), "SMV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "RC2"), "SMV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "BLF"), "RRV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "BNS"), "RRV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "ROS"), "RRV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "RAD"), "SRH", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "BOR"), "AV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "FAL"), "AV", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "GME"), "OR", ava)
  ava <- ifelse(str_detect(string = ava, pattern = "ZNW"), "OR", ava)

  
  ## separate years into train and validation sets
  ## semi-randomly create test set
  train_ind <- createDataPartition(ava, p = .7, list = FALSE, times = 1)
  train <- counts[train_ind, ]
  test <- counts[-train_ind, ]
  train$ava <- ava[train_ind]
  test$ava <- ava[-train_ind]
  
  ## tune/optimize model
  tmp <- train
  colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
  task <- makeClassifTask(data = tmp, target = "ava") # make an mlr task with counts17
  res <- tuneRanger(task, num.threads = 3)             # run tuning process
  write_tsv(res$recommended.pars, rec_pars_file)       # write model parameters to a file
  
  ## extract model parameters and use to build an optimal RF
  train$ava <- as.factor(train$ava)                # convert response to a factor
  optimal_rf <- ranger(
    dependent.variable.name = "ava",
    mtry            = res$recommended.pars$mtry,
    num.trees       = 10000,
    data            = train,
    sample.fraction = res$recommended.pars$sample.fraction,
    min.node.size   = res$recommended.pars$min.node.size,
    seed             = 1,
    importance       = 'permutation',
    local.importance = T
  )
  saveRDS(optimal_rf, file = optimal_rf_rds)
  
  # evaluate the accuracy of the model and generate a confusion matrix
  evaluate_model(optimal_ranger = optimal_rf,  data = train, reference_class = train$ava, 
                 plt_title = "Training Performance") 
  # validation data
  test$ava <- as.factor(test$ava)
  evaluate_model(optimal_ranger = optimal_rf, data = test, reference_class = test$ava, 
                 plt_title = "Validation Performance") 
  
  return(optimal_rf)
}

ava_train70 <- split_ava_model(counts = bac_dada,
                              rec_pars_file = "amplicon/bac/ava_train70_rec_pars.tsv",
                              optimal_rf_rds = "amplicon/bac/optimal_rf_ava_train70.RDS")

# predicted
# observed AS AV CRN OR RRV SMV SNC SRH
# AS   6  1   0  0   0   0   0   0
# AV   0  5   0  0   0   0   0   0
# CRN  0  0   3  0   0   0   0   0
# OR   0  0   0  6   0   0   0   0
# RRV  0  0   0  0   9   0   0   0
# SMV  0  0   0  0   0   7   0   0
# SNC  0  0   0  0   2   0   4   0
# SRH  0  0   0  0   0   0   0   3
# [1] "ACCURACY = 0.934782608695652"