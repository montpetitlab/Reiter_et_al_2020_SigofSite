library(ranger)
library(dplyr)
library(readr)

remotes::install_github("silkeszy/Pomona")
library(Pomona)
file.create(snakemake@output[['pomona']])