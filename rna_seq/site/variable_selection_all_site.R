library(ranger)
library(dplyr)
library(readr)
library(Pomona)
library(tibble)

set.seed(1)

# counts <- read_tsv("summarized_counts.tsv") # read in hash abund table
counts <- readRDS("rna_seq/site/summarized_counts.RDS") %>%
  as.data.frame(counts) %>%         # transform to dataframe
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
                         no.threads = 2, 
                         method = "ranger", type = "classification")
saveRDS(vita, "rna_seq/site/site_vita_all.RDS")

# write files -------------------------------------------------------------

## write predictive hashes
var <- vita$var                 # separate out selected predictive genes
write.table(var, "rna_seq/site/site_vita_all_vars.txt", quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
filt <- counts[ , colnames(counts) %in% var] # subset ibd to hashes in ibd_vita
write.csv(filt, "rna_seq/site/site_vita_all_filt.csv", quote = F)


# check output ----------------------------------------------------------------

# AWRI
# BCIN
# C5L36
# CUO98
# D499
# D6D18
# KLTH0C
# LPNK
# M438
# metfru2
# MWSF
# PEGC
# VIT
# Y/t/sn


#AWRI3578          BcDW1           BCIN          C5L36          CU098           D499          D6D18 
#186             18             71             36            354            419            522 
tmp <- gsub("_.*", "", var)
tmp <- gsub("KLTH0.*", "KLTH0", tmp)
tmp <- gsub("PEGC.*", "PEGC", tmp)
tmp <- gsub("MWSF.*", "MWSF", tmp)
tmp <- gsub("^Y.*", "scer", tmp)
tmp <- gsub("sn.*", "scer", tmp)
tmp <- gsub("LPNK.*", "LPNK", tmp)
tmp <- gsub("^t.*", "scer", tmp)
tmp <- gsub("RPR1", "scer", tmp)
tmp <- gsub("SCR1", "scer", tmp)
tmp <- gsub("ZOD1", "scer", tmp)
tmp <- gsub("RDN5.1", "scer", tmp)
tmp <- gsub("RUF22", "scer", tmp)
tmp <- gsub("HRA1", "scer", tmp)
tmp <- gsub("IRT1", "scer", tmp)
tmp <- gsub("TLC1", "scer", tmp)
tmp <- gsub("LSR1", "scer", tmp)

table(tmp)
