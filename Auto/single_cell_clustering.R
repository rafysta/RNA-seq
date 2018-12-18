#!/usr/bin/Rscript
# Automation of clustering of single cell RNA-seq 

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="file of normalized matrices"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-m", "--meta"), default="NA", help="meta file information for each cell. 1 clolumn should have Cell"),
  make_option(c("--color"), default="NA", help="custome colors. separated by ,"),
  make_option(c("--gene"), default="NA", help="file of gene list"),
  make_option(c("--cell"), default="NA", help="file of cell list")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mat <- as.character(opt["in"])
DIR_out <- as.character(opt["out"])
if(substring(DIR_out, nchar(DIR_out), nchar(DIR_out)) != "/"){
  DIR_out <- paste0(DIR_out, "/")
}
FILE_meta <- as.character(opt["meta"])

FILE_serat <- paste0(DIR_out, "serat.rds")


if (!requireNamespace("Rafysta", quietly = TRUE)){
  devtools::install_github("rafysta/rafysta")
}
library(Rafysta)

### Color palette
Colors <- as.character(opt["color"])
if(Colors != "NA"){
  Colors <- unlist(strsplit(Colors, ","))
}else{
  Colors <- NULL
}

FILE_gene <- as.character(opt["gene"])
if(FILE_gene != "NA"){
  Gene_list <- read.table(FILE_gene, header=FALSE, sep="\t", stringsAsFactors = FALSE)
  Gene_list <- Gene_list[,1]
}else{
  Gene_list <- NULL
}
FILE_cell <- as.character(opt["cell"])
if(FILE_cell != "NA"){
  Cell_list <- read.table(FILE_cell, header=FALSE, sep="\t", stringsAsFactors = FALSE)
  Cell_list <- Cell_list[,1]
}else{
  Cell_list <- NULL
}

#=============================================
# Scaling by Seurat
#=============================================
suppressWarnings(suppressMessages(library(Seurat)))
if(!file.exists(FILE_serat)){
  mat <- readRDS(FILE_mat)
  
  # Seratで解析。normalizeはしない
  serat <- CreateSeuratObject(raw.data = mat)
  rm(mat)
  
  # scale
  serat <- ScaleData(object = serat, genes.use = Gene_list,  vars.to.regress = c("nUMI"))
  
  saveRDS(serat, FILE_serat)
}else{
  serat <- readRDS(FILE_serat)
}


# Cell information
if(FILE_meta != "NA"){
  D_cell <- read.table(FILE_meta, header=TRUE, sep="\t", stringsAsFactors = FALSE)
}


#=============================================
# PCA
#=============================================
FILE_pca <- paste0(DIR_out, "pca.rds")
if(!file.exists(FILE_pca)){
  pca <- Clustering_PCA(serat, genes=Gene_list, cells=Cell_list)
  saveRDS(pca, FILE_pca)
}else{
  pca <- readRDS(FILE_pca)
}
plot_PCA(pca, file=paste0(DIR_out, "pca.png"))


ok_cell <- Filter_PCA(pca, file=paste0(DIR_out, "pca_outlier.png"))
FILE_pca2 <- paste0(DIR_out, "pca2.rds")
if(!file.exists(FILE_pca2)){
  pca2 <- Clustering_PCA(serat, genes=Gene_list, cells=ok_cell)
  saveRDS(pca2, FILE_pca2)
}else{
  pca2 <- readRDS(FILE_pca2)
}
p1 <- plot_PCA(pca2, file=paste0(DIR_out, "pca2.png"))



#=============================================
# tSNE
#=============================================
FILE_tsne <- paste0(DIR_out, "tsne.rds")
if(!file.exists(FILE_tsne)){
  tsne <- Clustering_tSNE(pca, pca_use = 1:20)
  saveRDS(tsne, FILE_tsne)
}else{
  tsne <- readRDS(FILE_tsne)
}
p2 <- plot_tSNE(tsne, file=paste0(DIR_out, "tsne.png"))



#=============================================
# SIMILR
#=============================================
FILE_SIMILR <- paste0(DIR_out, "sim.rds")
if(!file.exists(FILE_SIMILR)){
  sim <- Clustering_SIMILR(pca)
  saveRDS(sim, FILE_SIMILR)
}else{
  sim <- readRDS(FILE_SIMILR)
}
p3 <- plot_SIMILR(sim, file=paste0(DIR_out, "SIMILR.png"))


#=============================================
# Mix three graph
#=============================================
p <- plot_grid(p1, p2, p3, nrow=3, align="v")
save_plot(paste0(DIR_out, "clustering_mix.png"), p, ncol=1, nrow=3)



