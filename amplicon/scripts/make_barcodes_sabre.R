# make adapter table for demultiplexing with sabre
barcodes <- read.csv(snakemake@input[['primers']], stringsAsFactors = F)
barcodes <- barcodes[ , c("sequence_name", "blendID")]
colnames(barcodes) <- c("barcode", "r1")
barcodes$barcode <- gsub("CTACCTGCGGARGGATCA", "", barcodes$barcode)
barcodes$r1 <- gsub("Envir. Neg ", "negctrl", barcodes$r1)
barcodes$r1 <- gsub("mock community #", "mockcomm", barcodes$r1)
barcodes$r1 <- gsub("negative control #", "negctrl", barcodes$r1)
barcodes$r1 <- gsub("negative control#", "negctrl", barcodes$r1)
barcodes$r2 <- paste0(barcodes$r1, "_R2.fq")
barcodes$r1 <- paste0(barcodes$r1, "_R1.fq")

write.table(barcodes, snakemake@output[['barcodes']], 
            row.names = F, quote = F, sep = "\t", col.names = F)