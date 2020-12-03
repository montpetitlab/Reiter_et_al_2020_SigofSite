library(ranger)
library(dplyr)
library(readr)
library(Pomona)
library(tibble)

set.seed(1)
# counts <- summarized_counts
counts <- readRDS("rna_seq/avas/summarized_counts_ava.RDS") %>% # read in counts
  as.data.frame() %>%            # transform to dataframe
  column_to_rownames("sample")   # set sample as rownames and rm

# head(colnames(counts))
## make classification vector
ava <- counts$ava

# subtract AVA
counts <- counts[ , -1]


# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
vita <- var.sel.vita(x = counts, y = ava, p.t = 0.05,
                     ntree = 5000, mtry.prop = 0.2, nodesize.prop = 0.1,
                     no.threads = 3, 
                     method = "ranger", type = "classification")
saveRDS(vita, "rna_seq/avas/ava_vita_all.RDS")

# write files -------------------------------------------------------------

## write predictive hashes
var <- vita$var                 # separate out selected predictive genes
write.table(var, "rna_seq/avas/ava_vita_all_vars.txt", quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
filt <- counts[ , colnames(counts) %in% var] # subset ibd to hashes in ibd_vita
write.csv(filt, "rna_seq/avas/ava_vita_all_filt.csv", quote = F)