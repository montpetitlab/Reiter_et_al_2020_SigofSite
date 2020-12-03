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

seed <- as.numeric(snakemake@wildcards[['seed']])
set.seed(seed)

#filt <- read_csv("rna_seq/ava_seeds/ava_vita_all_filt_seed1.csv") %>%
filt <- read_csv(snakemake@input[['counts_filt']]) %>%
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

# tune ranger -------------------------------------------------------------

# Make an mlr task with the ibd_train dataset here 
tmp <- train
colnames(tmp) <-  make.names(colnames(tmp))
task <- makeClassifTask(data = tmp, target = "ava")
# Run tuning process
res <- tuneRanger(task, num.threads = 3)

# write model parameters to a file
write_tsv(res$recommended.pars, snakemake@output[['rec_pars']])

# build optimal model ----------------------------------------------------------

# extract model parameters and use to build an optimal RF

# use model parameters to build optimized RF
train$ava <- as.factor(train$ava)
optimal_rf <- ranger(
  dependent.variable.name = "ava",
  mtry            = res$recommended.pars$mtry,
  num.trees       = 10000,
  data            = train,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed            = seed,
  importance       = 'permutation',
  local.importance = T
)

saveRDS(optimal_rf, file = snakemake@output[['optimal_rf']])

plt <- evaluate_model(optimal_ranger = optimal_rf,
                       data = test, reference_class = as.factor(test$ava),
                       plt_title = "Test Performance")

ggsave(filename = snakemake@output[['confusion_plt']], plot = plt, 
       scale = 1, width = 6, height = 4, dpi = 300)