# Evaluates RF performance on common set of genes attained by running variable
# selection on all samples.

library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
library(tibble)
source("scripts/evaluate_model.R")
source("scripts/ggplotConfusionMatrix.R")

set.seed(1)

filt <- read_csv("rna_seq/avas/ava_vita_all_filt.csv") %>%
  as.data.frame() %>%
  column_to_rownames("X1")

## separate years
filt19 <- grepl(pattern = "2019", x = rownames(filt))
filt19 <- filt[filt19, ]
filt17 <- grepl(pattern = "2017", x = rownames(filt))
filt17 <- filt[filt17, ]

filt19$ava <- gsub("_.*", "", rownames(filt19))
filt17$ava <- gsub("_.*", "", rownames(filt17))


# train 2017, validate on 2019 --------------------------------------------

# tune
tmp <- filt17
colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
task <- makeClassifTask(data = tmp, target = "ava")  # make an mlr task with filt17
res <- tuneRanger(task, num.threads = 3)             # run tuning process
write_tsv(res$recommended.pars, "rna_seq/avas/ava_train2017_rec_pars.tsv") # write model parameters to a file

# extract model parameters and use to build an optimal RF
filt17$ava <- as.factor(filt17$ava)                # convert response to a factor
optimal_rf17 <- ranger(
  dependent.variable.name = "ava",
  mtry            = res$recommended.pars$mtry,
  num.trees        = 10000,
  data             = filt17,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed             = 2,
  importance       = 'permutation',
  local.importance = T
)

# View(optimal_rf17$variable.importance.local)
saveRDS(optimal_rf17, file = "rna_seq/avas/ava_train2017_optimal_rf.RDS")

# evaluate the accuracy of the model and generate a confusion matrix
evaluate_model(optimal_ranger = optimal_rf17,  data = filt17, reference_class = filt17$ava, 
               plt_title = "2017 Training Performance") # eval training data -- 2017
# validation data
evaluate_model(optimal_ranger = optimal_rf17, data = filt19, reference_class = as.factor(filt19$ava),
               plt_title = "2019 Validation Performance") # eval 2017 model on 2019 data

# 57% accuracy on 2019 data

# train 2019, validate on 2017 --------------------------------------------

# tune
tmp <- filt19
colnames(tmp) <-  make.names(colnames(tmp))          # tmp change names to make compatible with tuning
task <- makeClassifTask(data = tmp, target = "ava")  # make an mlr task with filt19
res <- tuneRanger(task, num.threads = 3)             # run tuning process
write_tsv(res$recommended.pars, "rna_seq/avas/ava_train2019_rec_pars.tsv") # write model parameters to a file

# extract model parameters and use to build an optimal RF
filt19$ava <- as.factor(filt19$ava)                # convert response to a factor
optimal_rf19 <- ranger(
  dependent.variable.name = "ava",
  mtry             = res$recommended.pars$mtry,
  num.trees        = 10000,
  data             = filt19,
  sample.fraction  = res$recommended.pars$sample.fraction,
  min.node.size    = res$recommended.pars$min.node.size,
  seed             = 1,
  importance       = 'permutation',
  local.importance = T
)

saveRDS(optimal_rf19, file = "rna_seq/avas/ava_train2019_optimal_rf.RDS")

# evaluate the accuracy of the model and generate a confusion matrix
evaluate_model(optimal_ranger = optimal_rf19,  data = filt19, reference_class = filt19$ava, 
               plt_title = "2019 Training Performance") # eval training data -- 2019
# validation data
evaluate_model(optimal_ranger = optimal_rf19, data = filt17, reference_class = as.factor(filt17$ava), 
               plt_title = "2017 Validation Performance") # eval 2019 model on 2017 data

# 63.3% accuracy
# eval variable importance/local importance -------------------------------

which.max(optimal_rf17$variable.importance.local[1, ])
which.max(optimal_rf17$variable.importance.local[26, ])


var17 <- data.frame(gene = names(optimal_rf17$variable.importance), imp = optimal_rf17$variable.importance)
var17 <- var17 %>%
  filter(imp >0) %>%
  arrange(desc(imp))

var19 <- data.frame(gene = names(optimal_rf19$variable.importance), imp = optimal_rf19$variable.importance)
var19 <- var19 %>%
  filter(imp >0) %>%
  arrange(desc(imp))

table(var19$gene[1:100] %in% var17$gene[1:100])
table(var19$gene %in% var17$gene)

