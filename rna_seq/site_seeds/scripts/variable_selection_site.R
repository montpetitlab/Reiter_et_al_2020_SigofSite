library(ranger)
library(dplyr)
library(readr)
library(Pomona)
library(tibble)

set.seed(as.numeric(snakemake@wildcards[['seed']]))

counts <- readRDS(snakemake@input[['counts']]) %>%
  as.data.frame() %>%               # transform to dataframe
  column_to_rownames("ts_sum_id")   # set sample as rownames

## make classification vector
site <- gsub("201[679]{1}_", "", rownames(counts))
site <- gsub("_.*", "", site)


# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
vita <- var.sel.vita(x = counts, y = site, p.t = 0.05,
                         ntree = 5000, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = 3, 
                         method = "ranger", type = "classification")
saveRDS(vita, snakemake@output[['vita_rds']])

# write files -------------------------------------------------------------

## write predictive hashes
var <- vita$var                 # separate out selected predictive genes
write.table(var, snakemake@output[['vita_vars']], quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
filt <- counts[ , colnames(counts) %in% var] # subset ibd to hashes in vita
write.csv(filt, snakemake@output[["counts_filt"]], quote = F)