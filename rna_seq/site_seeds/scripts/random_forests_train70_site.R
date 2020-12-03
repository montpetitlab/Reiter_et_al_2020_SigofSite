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

filt <- read_csv(snakemake@input[['counts_filt']]) %>%
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

# Make an mlr task with the train dataset here 
tmp <- train
colnames(tmp) <-  make.names(colnames(tmp))
task <- makeClassifTask(data = tmp, target = "site")
# Run tuning process
res <- tuneRanger(task, num.threads = 3)

# write model parameters to a file
write_tsv(res$recommended.pars, snakemake@output[['rec_pars']])

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
  seed             = seed,
  importance       = 'permutation',
  local.importance = T,
  keep.inbag       = T
)


saveRDS(optimal_rf, file = snakemake@output[['optimal_rf']])

plt <- evaluate_model(optimal_ranger = optimal_rf,
                      data = test, reference_class = as.factor(test$site),
                      plt_title = "Test Performance")

ggsave(filename = snakemake@output[['confusion_plt']], plot = plt,
       scale = 1, width = 6, height = 4, dpi = 300)
