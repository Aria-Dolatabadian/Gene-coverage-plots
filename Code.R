#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")



#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome")



library(GenVisR)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm9)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(BSgenome)


# read the coverage data into R
covData <- read.delim("STAT1_mm9_coverage.tsv")


#Formating coverage data

# rename the columns
colnames(covData) <- c("chromosome", "start", "end", "TAC245", "TAC265")

# create a function to split the data frame into lists of data frames
samples <- c("TAC245", "TAC265")
a <- function(x, y){
    col_names <- c("chromosome", "end", x)
    y <- y[,col_names]
    colnames(y) <- c("chromosome", "end", "cov")
    return(y)
}
covData <- lapply(samples, a, covData)

names(covData) <- samples


library("BSgenome")
available.genomes()


library("BSgenome.Mmusculus.UCSC.mm9")
genomeObject <- BSgenome.Mmusculus.UCSC.mm9

library("TxDb.Mmusculus.UCSC.mm9.knownGene")
TxDbObject <- TxDb.Mmusculus.UCSC.mm9.knownGene

#Creating a Granges object

# get the chromosome, and the start and end of our coverage data
chromosome <- as.character(unique(covData[[1]]$chromosome))
start <- as.numeric(min(covData[[1]]$end))
end <- as.numeric(max(covData[[1]]$end))

# define the genomic range
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start, end=end))

#Creating an initial coverage plot

# create an initial plot
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line")



# add a flank and create another plot
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start-500, end=end+500))
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line", gene_isoformSel="uc007aya.1")


# adjust compression of genomic features
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line", gene_isoformSel="uc007aya.1", base=NA, transform=NA)







