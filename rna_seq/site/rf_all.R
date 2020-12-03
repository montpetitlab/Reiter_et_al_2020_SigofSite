library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
library(caret)
library(tibble)
# source(snakemake@input[['eval_model']])
# source(snakemake@input[['ggconfusion']])

set.seed(1)

filt <- read_csv("rna_seq/site/site_vita_all_filt.csv") %>%
  as.data.frame() %>%
  column_to_rownames("X1")

## semi-randomly create test set
site <- gsub("201[79]{1}_", "", rownames(filt))
site <- gsub("_.*", "", site)

train_ind<- createDataPartition(site, 
                                p = .7,
                                list = FALSE,
                                times = 1)
train <- filt[train_ind, ]
test <- filt[-train_ind, ]

# remove testdata from training data
train$site <- gsub("201[79]{1}_", "", rownames(train))
train$site <- gsub("_.*", "", train$site)
test$site <- gsub("201[79]{1}_", "", rownames(test))
test$site <- gsub("_.*", "", test$site)

# check that all sites are represented in test data
table(train$site)
table(test$site)
# tune ranger -------------------------------------------------------------

# Make an mlr task with the ibd_train dataset here 
tmp <- train
colnames(tmp) <-  make.names(colnames(tmp))
task <- makeClassifTask(data = tmp, target = "site")
# Run tuning process
res <- tuneRanger(task, num.threads = 3)

# write model parameters to a file
write_tsv(res$recommended.pars, "rna_seq/site/site_train70_rec_pars.tsv")

# build optimal model ----------------------------------------------------------

# extract model parameters and use to build an optimal RF

# use model parameters to build optimized RF
train$site <- as.factor(train$site)
optimal_rf <- ranger(
  dependent.variable.name = "site",
  mtry             = res$recommended.pars$mtry,
  num.trees        = 10000,
  data             = train,
  sample.fraction  = res$recommended.pars$sample.fraction,
  min.node.size    = res$recommended.pars$min.node.size,
  seed             = 5,
  importance       = 'permutation',
  local.importance = T,
  keep.inbag       = T
)


saveRDS(optimal_rf, file = 'rna_seq/site/site_train70_optimal_rf.RDS')
# evaluate the accuracy of the model and generate a confusion matrix
# training data
evaluate_model(optimal_ranger = readRDS("rna_seq/site/site_train70_optimal_rf.RDS"),  
               data = train, reference_class = as.factor(train$site), 
               plt_title = "Training Performance")

evaluate_model(optimal_ranger = readRDS("rna_seq/site/site_train70_optimal_rf.RDS"),  
               data = test, reference_class = as.factor(test$site), 
               plt_title = "Test Performance")
