# counts matrix and DESeq2 on TARGET counts matrix


# assign path to folder in which I have permission to do stuff

myPaths <- .libPaths()   # get the paths
.libPaths(myPaths)  # reassign them

# -------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

if ("dplyr" %in% rownames(install.packages())==FALSE)
  install.packages('dplyr')





library(dplyr)
library(DESeq2)



files <- list.files('../data/counts/')

# option 1: cbind --------------

# create first dataframe


filename = paste0('../data/counts/',files[1])

dataset <- read.table(filename,header = FALSE, row.names = 1)
dataset <- dataset[5:nrow(dataset),]
dataset <- dataset['V2']
colnames(dataset) <- sub('_.*','',files[1])



for (f in files[2:length(files)]){
  
  # Create the first data if no data exist yet
  
  filename = paste0('../data/counts/',f)
  
  # if data already exist, then append it together
  
  tempory <- read.table(filename, header=FALSE,row.names = 1)
  tempory <- tempory[5:nrow(tempory),]
  tempory <- tempory['V2']
  colnames(tempory) <- sub('_.*','',f) # delete ReadsPerGene ending
  dataset <- cbind(dataset,tempory)
  rm(tempory)
  
}

dataset.sorted <- dataset[order(row.names(dataset)),]


write.csv(dataset.sorted,file='../data/tcga_unnormalized_counts.csv')



# normalize counts matrix with deseq2 -----

cts <- as.matrix(dataset.sorted)

coldata <- c(condition = 'Standard') 
coldata <- rep(coldata,136)

coldata <- data.frame(coldata)


coldata$coldata <- as.factor(coldata$coldata)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 1)


dds <- DESeq(dds)

res <- results(dds)


res.frame <- as.data.frame(res)

counts_normalized = counts(dds,normalized=TRUE);
counts_normalized.log <- log2(1+counts_normalized)



sizeFactors <- sizeFactors(dds)

write.csv(sizeFactors,"../supplementary/sizeFactors.csv")


# remove rows with all zeros ---------

counts_normalized <- counts_normalized[rowSums(counts_normalized) > 0, ]
counts_normalized.log <- counts_normalized.log[rowSums(counts_normalized.log) > 0, ]



write.csv(counts_normalized,'../data/tcga_normalized_counts.csv')
write.csv(counts_normalized.log, '../data/tcga_log_normalized_counts.csv')





