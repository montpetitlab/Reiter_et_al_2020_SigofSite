library(phylotools)

# make adapter fasta for demultiplexing with cutadapt
barcodes <- read.csv(snakemake@input[['primers']], stringsAsFactors = F)
barcodes <- barcodes[ , c("blendID", "sequence_name")]
colnames(barcodes) <- c("seq.name", "seq.text")
barcodes$seq.name <- gsub("Envir. Neg ", "negctrl", barcodes$seq.name)
barcodes$seq.name <- gsub("mock community #", "mockcomm", barcodes$seq.name)
barcodes$seq.name <- gsub("negative control #", "negctrl", barcodes$seq.name)
barcodes$seq.name <- gsub("negative control#", "negctrl", barcodes$seq.name)
dat2fasta(barcodes, outfile = snakemake@output[['barcodes']])