# Longitudinal single cell RNA-seq of the inflamed colon
# Mice were treated with 1.5%DSS for 6 days then switched to drinking water. 
# Colon from each time points (Day 0, 3, 6, 9, 12 and 15) or treatment 
# (PBS or Serpina3n) were collected. To limit mouse-to-mouse biased, 
# three biological replicates from each group were pooled and analyzed 
# using FACS-based, smart-seq2 method. The sequencing libraries were
# sequenced on the Illumina NextSeq500 platform.
library(Seurat)
library(stringr)
library(data.table)

setwd("~/data_Albas/Mice_Paper")

counts <- as.data.frame(fread('GSE148794_tc_ibd.count_table.tsv', sep = '\t'))
rownames(counts) <- counts$V1
counts <- counts[,-grep("V1", colnames(counts))]
counts <- as.matrix(counts)

metadata <- read.delim('GSE148794_tc_ibd.metadata.tsv', sep = '\t')
ref <- gsub('\\.','_',paste0('R19', str_split_fixed(pattern = '.19', str = metadata$Sample_name, n = 2)[,2]))
refok <- str_split_fixed(string = ref, pattern = '_S\\d', n = 2)[,1]

timepoint <- gsub(pattern = 'tc_ibd.',replacement = '',
                  x = str_split_fixed(pattern = '.19', 
                                      str = metadata$Sample_name,
                                      n = 2)[,1])

metadata$refok <- refok
metadata$timepoint <- timepoint

dat <- data.frame()
for(colname in colnames(counts)){
  ref <- substr(colname, 0,42)
  nn <-grep(ref, metadata$refok)
  dat <- rbind(dat, metadata[nn,])
  rownames(dat)[nrow(dat)] <- colname
}

saveRDS(dat, file = 'GSE148794_metadata_ok.RDS')


