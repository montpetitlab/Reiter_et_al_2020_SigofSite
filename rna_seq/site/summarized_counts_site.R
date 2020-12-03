library(dplyr)
library(readr)
library(tidyr)

setwd("~/github/2020-pn-site")


# read in samples ---------------------------------------------------------
## info
info <- read_csv("samples_ava.csv") %>%
  filter(!is.na(tank)) %>%
  filter(hr >10) %>%
  filter(hr < 38 | hr > 46) %>%
  filter(hr <74 | hr >89) %>%
  select(id, ts_sum_id, year, tank, ava_id)

# re-add MSA and P5A
info2 <- read_csv("samples_ava.csv") %>%
  filter(ava_id %in% c("AS1", "AS2")) %>%
  filter(year == 2017) %>%
  filter(hr > 70) %>%
  select(id, ts_sum_id, year, tank, ava_id)

info <- rbind(info2, info)

## 2019
v19 <- read_tsv("rna_seq/v2019/outputs/counts_all_organisms/raw_counts.tsv") %>%
  filter(!gene %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"))
colnames(v19) <- gsub("outputs\\/htseq_all_organisms\\/", "", colnames(v19))
colnames(v19) <- gsub("_readcounts\\.txt", "", colnames(v19))
v19 <- select(v19, c(colnames(v19)[colnames(v19) %in% info$id], gene))

v17 <- read_tsv("rna_seq/v2017/outputs/counts_all_organisms/raw_counts.tsv") %>%
  filter(!gene %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"))
colnames(v17) <- gsub("outputs\\/htseq_all_organisms\\/", "", colnames(v17))
colnames(v17) <- gsub("_readcounts\\.txt", "", colnames(v17))
v17 <- select(v17, c(colnames(v17)[colnames(v17) %in% info$id], gene))

# join on gene, remove rows that sum to zero, normalize by library size, wide -> long, join with info, group_by year/tank/site/gene, 
# calculate mean, min, max, 

counts <- full_join(v17, v19, by = "gene")
counts <- counts[rowSums(counts[-which(names(counts) %in% "gene")]) > 0, ] # remove rows that sum to zero
counts <- counts %>%
  mutate_at(vars(-gene), funs(./sum(.)))
counts <- pivot_longer(counts, cols = -gene, names_to = 'id', values_to = 'count')
counts <- left_join(counts, info, by = "id")

summarized_counts <- counts %>%
  group_by(ts_sum_id, gene) %>%
  summarize(mean = mean(count, na.rm = TRUE),
            min = min(count),
            max = max(count),
            raw = sum(count),
            deviation = sd(count))

summarized_counts <- summarized_counts %>%
  ungroup() %>%
  pivot_longer(cols = mean:deviation, names_to = "metric", values_to = "value") %>%
  mutate(gene = paste(gene, metric, sep = "_")) %>%
  select(-metric)

summarized_counts <- pivot_wider(summarized_counts, id_cols = ts_sum_id, names_from = gene, values_from = value)

saveRDS(summarized_counts, "rna_seq/site/summarized_counts2.RDS")
write_tsv(summarized_counts, "rna_seq/site/summarized_counts2.tsv")
