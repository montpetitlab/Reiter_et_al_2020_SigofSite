library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
library(caret)
library(tibble)
source("scripts/evaluate_model.R")
source("scripts/ggplotConfusionMatrix.R")

set.seed(3)

filt <- read_csv("rna_seq/avas/ava_vita_all_filt.csv") %>%
  as.data.frame() %>%
  column_to_rownames("X1")

## semi-randomly create test set
train_ind<- createDataPartition(gsub("_.*", "", rownames(filt)), 
                                p = .7,
                                list = FALSE,
                                times = 1)

train <- filt[train_ind, ]
test <- filt[-train_ind, ]

# remove testdata from training data
train$ava <- gsub("_.*", "", rownames(train))
test$ava <- gsub("_.*", "", rownames(test))

# check that all sites are represented in test data
table(train$ava)
table(test$ava)
# tune ranger -------------------------------------------------------------

# Make an mlr task with the ibd_train dataset here 
tmp <- train
colnames(tmp) <-  make.names(colnames(tmp))
task <- makeClassifTask(data = tmp, target = "ava")
# Run tuning process
res <- tuneRanger(task, num.threads = 3)

# write model parameters to a file
write_tsv(res$recommended.pars, "rna_seq/avas/ava_train70_rec_pars.tsv")

# build optimal model ----------------------------------------------------------

# extract model parameters and use to build an optimal RF

# use model parameters to build optimized RF
train$ava <- as.factor(train$ava)
optimal_rf <- ranger(
  dependent.variable.name = "ava",
  mtry            = res$recommended.pars$mtry,
  num.trees       = 15000,
  data            = train,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed            = 1,
  importance       = 'permutation',
  local.importance = T
)

saveRDS(optimal_rf, file = 'rna_seq/avas/ava_train70_optimal_rf.RDS')

evaluate the accuracy of the model and generate a confusion matrix
training data
evaluate_model(optimal_ranger = readRDS("rna_seq/avas/ava_train70_optimal_rf.RDS"),
               data = train, reference_class = train$ava,
               plt_title = "Training Performance")

evaluate_model(optimal_ranger = readRDS("rna_seq/avas/ava_train70_optimal_rf.RDS"),
               data = test, reference_class = as.factor(test$ava),
               plt_title = "Test Performance")

# eval --------------------------------------------------------------------

View(optimal_rf$variable.importance)

