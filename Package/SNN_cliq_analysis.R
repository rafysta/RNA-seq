# SNN-Cliq analysis

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="read count object"),
  make_option(c("-L", "--minimum_read"), default=118000, help="minimum required reads"),
  make_option(c("-o", "--out"),help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_object <- as.character(opt["in"])
FILE_OUT <- as.character(opt["out"])
Threshold_minRead <- as.numeric(as.character(opt["minimum_read"]))


DATA <- readRDS(FILE_object)
nameList <- strsplit(colnames(DATA), ":")
nameMatrix <- matrix(unlist(nameList), ncol=2, byrow=TRUE)
# table(as.character(nameMatrix[,1]))


# Normalize read number for each cells.
TotalRead_per_cells <- apply(DATA, 2, sum)
AverageRead_per_cells <- mean(TotalRead_per_cells)
Scale_factor <- AverageRead_per_cells / TotalRead_per_cells
DATA.normalized <- t(t(DATA) * Scale_factor)
DATA.normalized.log10 <- log10(DATA.normalized + 1)

# read threshold
cell_target <- which(TotalRead_per_cells > Threshold_minRead)


dataForPCA <- t(DATA.normalized.log10[,cell_target])


#### SNN_Cliq
source("/applications/SNNCliq/current/SNN.R")
SNN(dataForPCA, FILE_OUT, k=3, distance="euclidean")











