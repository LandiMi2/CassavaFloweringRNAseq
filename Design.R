library(data.table)
library(DESeq2)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Flowering/counts/")

counts <- fread("counts.txt")

#remove annotation info
counts <- counts[, -c("Chr", "Start", "End", "Strand", "Length")]

#now rename the sample names 
colnames(counts) <- sub(
  "^.*\\/(T[0-9]+_[0-9]+)_Aligned.*$",
  "\\1",
  colnames(counts)
)

##filtering step - for all (let's see)
# keeps genes that have a count of >10 in at least 6 samples
sample.cols <- setdiff(names(counts), "Geneid")
#filter
keep <- rowSums(counts[, ..sample.cols] >= 10) >= 6
counts <- counts[keep]
counts

###count matrix
count.mat <- as.matrix(counts[, ..sample.cols])
rownames(count.mat) <- counts$Geneid

#colDARA
sample.names <- colnames(count.mat)

colData <- data.frame(
  row.names = sample.names,
  time = sub("^(T[0-9]+)_.*", "\\1", sample.names),
  rep  = as.integer(sub(".*_(\\d+)$", "\\1", sample.names))
)

#light vs no lit
colData$condition <- ifelse(colData$rep <= 9, "Light", "NoLight")

#add genotype 
colData$geno <- NA #empty colum

#lit samples (1–9)
colData$geno[colData$rep %in% 1:3] <- "DSC120_late"
colData$geno[colData$rep %in% 4:6] <- "DSC196_non"
colData$geno[colData$rep %in% 7:9] <- "DSC272_early"

# NoLight samples (10–18)
colData$geno[colData$rep %in% 10:12] <- "DSC120_late"
colData$geno[colData$rep %in% 13:15] <- "DSC196_non"
colData$geno[colData$rep %in% 16:18] <- "DSC272_early"


#order -as factor -levels
colData$time <- factor(colData$time, levels = c("T0","T1","T2"))
colData$condition <- factor(colData$condition, levels = c("NoLight","Light"))
colData$geno <- factor(colData$geno,
                       levels = c("DSC120_late","DSC196_non","DSC272_early"))

#relevel - I want the early genotype to be the reference
colData$geno <- relevel(colData$geno, ref = "DSC272_early")
colData$time <- relevel(colData$time, ref = "T0")
colData$condition <- relevel(colData$condition, ref = "NoLight")

colData

#all sample present
table(colData$condition, colData$geno, colData$time)

#Now lets run deseq
dds <- DESeqDataSetFromMatrix(
  countData = count.mat,
  colData   = colData,
  design    = ~ condition * geno * time
)

## run DEGs
dds <- DESeq(dds)

##save results 
saveRDS(dds, file = "DESeqCasFlower.rds")
