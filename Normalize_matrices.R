#!/usr/bin/Rscript
# Normalize single cell RNA-seq matrix

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="input Robject(use FPKM table for TPM normalization"),
  make_option(c("--format"), default="counts", help="input file format is counts or fpkm"),
  make_option(c("-d", "--db"), default="NA", help="SQLite3 database"),
  make_option(c("-m", "--method"), default="scran", help="normalizing method, log, quantile, scran, cpm, upperquantile, tmm, TPM, center "),
  make_option(c("--cell"), default="NA", help="SQLite3 query to output Barcode list or file name of cell list(first column is cell)"),
  make_option(c("--gene"), default="NA", help="SQLite3 query to output gene list or file name of gene list(first column is gene)"),
  make_option(c("-o", "--out"), default="NA", help="output Robject"),
  make_option(c("--text"), default="NA", help="output text file"),
  make_option(c("-n", "--min_size"), default=30, help="minimum cell number for clustering for quick cluster by scran")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mat <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])

mat <- readRDS(FILE_mat)

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


suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
mat <- as.matrix(mat)


Down_Sample_Matrix <- function (expr_mat) {
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


switch(as.character(opt["method"]), 
   log={
     mat.norm <- log2(mat + 1) 
   },
   quantile={
     suppressPackageStartupMessages(library(preprocessCore))
     mat.norm <- normalize.quantiles(log2(mat + 1))
     rownames(mat.norm) <- rownames(mat)
     colnames(mat.norm) <- colnames(mat)
   },
   rle={
     sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
     rm(mat)
     mat.norm <- logcounts(normaliseExprs(
       sce,
       method = "RLE", 
       return_log = TRUE,
       return_norm_as_exprs = TRUE
     ))
   },
   scran={
     sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
     rm(mat)
     NUM_MIN_SIZE <- as.numeric(as.character(opt["min_size"]))
     clusters <- quickCluster(sce, min.size=NUM_MIN_SIZE)
     sce <- computeSumFactors(sce, cluster=clusters, positive=TRUE)
     sce <- normalize(sce)
     mat.norm <- logcounts(sce)
   },
   cpm={
     sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
     rm(mat)
     mat.norm <- log2(calculateCPM(sce, use.size.factors = FALSE) + 1) 
   },
  tmm={
    sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
    rm(mat)
    mat.norm <- logcounts(normaliseExprs(
      sce,
      method = "TMM",
      return_log = TRUE,
      return_norm_as_exprs = TRUE
    ))
   },
  tpm={
    if(as.character(opt["format"]) == "fpkm"){
      # follow Tirosh 2016 et al method. (divided by 10)
      mat.norm <- apply(mat, 2, fpkmToTpm)
      mat.norm <- log2(mat.norm/10 + 1)
    }else{
      sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
      rm(mat)
      mat.norm <- calculateTPM(
        sce,
        effective_length = 5e04,
        calc_from = "counts"
      )
      mat.norm <- log2(mat.norm/10 + 1)
    }
  },
  upperquantile={
    sce <- SingleCellExperiment(list(counts=as.matrix(mat)))
    rm(mat)
    mat.norm  <- logcounts(normaliseExprs(
      sce,
      method = "upperquartile", 
      p = 0.99,
      return_log = TRUE,
      return_norm_as_exprs = TRUE
    ))
  },
  downsample={
    mat.norm <- log2(Down_Sample_Matrix(mat) + 1)
  },
  center={
    # only centerizing (Tirosh 2016 et al) use normalized value
    mat.norm <- t(apply(mat, 1, scale, center=TRUE, scale=FALSE))
    colnames(mat.norm) <- colnames(mat)
  },
  {
    print('select correct normalization method')
    q()
  }
)


saveRDS(mat.norm, FILE_out)
FILE_text <- as.character(opt["text"])
if(FILE_text != "NA"){
  write.table(mat.norm, FILE_text, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
}




