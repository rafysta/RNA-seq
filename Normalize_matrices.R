#!/usr/bin/Rscript
# Normalize single cell RNA-seq matrix

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="input Robject(use FPKM table for TPM normalization"),
  make_option(c("-d", "--db"), default="NA", help="SQLite3 database"),
  make_option(c("-m", "--method"), default="scran", help="normalizing method, log quantile rle scran cpm upperquantile downsample tpm center "),
  make_option(c("--cell"), default="NA", help="SQLite3 query to output Barcode list or file name of cell list(first column is cell)"),
  make_option(c("--gene"), default="NA", help="SQLite3 query to output gene list or file name of gene list(first column is gene)"),
  make_option(c("-o", "--out"), default="NA", help="output R object"),
  make_option(c("--text"), default="NA", help="output text file"),
  make_option(c("-n", "--min_size"), default=30, help="minimum cell number for clustering for quick cluster by scran")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mat <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
METHOD <- as.character(opt["method"])
SIZE <- as.numeric(as.character(opt["min_size"]))

mat <- readRDS(FILE_mat)

if (!requireNamespace("Rafysta", quietly = TRUE)){
  devtools::install_github("rafysta/rafysta")
}
library(Rafysta)


### Cell list
if(as.character(opt["cell"]) != "NA"){
  if(file.exists(as.character(opt["cell"]))){
    D <- read.table(as.character(opt["cell"]), header=FALSE, sep="\t", stringsAsFactors = FALSE)
    CELLs <- D[,1]
  }else{
    library("RSQLite")
    FILE_DB <- as.character(opt["db"])
    DB <- dbConnect(SQLite(), FILE_DB)
    CELLs <- dbGetQuery(DB, as.character(opt["cell"]))[,1]
    dbDisconnect(DB)
  }
  CELLs <- intersect(CELLs, colnames(mat))
  mat <- mat[,CELLs]
}
### Gene list
if(as.character(opt["gene"]) != "NA"){
  if(file.exists(as.character(opt["gene"]))){
    D <- read.table(as.character(opt["gene"]), header=FALSE, sep="\t", stringsAsFactors = FALSE)
    GENEs <- D[,1]
  }else{
    library("RSQLite")
    FILE_DB <- as.character(opt["db"])
    DB <- dbConnect(SQLite(), FILE_DB)
    GENEs <- dbGetQuery(DB, as.character(opt["gene"]))[,1]
    dbDisconnect(DB)
  }
  GENEs <- intersect(GENEs, rownames(mat))
  mat <- mat[GENEs,]
}


### Normalize
mat.norm <- NormalizeMatrix(mat, method=METHOD, min_size=SIZE)


saveRDS(mat.norm, FILE_out)
FILE_text <- as.character(opt["text"])
if(FILE_text != "NA"){
  write.table(mat.norm, FILE_text, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
}




