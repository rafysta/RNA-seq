#!/usr/bin/Rscript
# Automation of clustering of single cell RNA-seq 

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="file of normalized matrices"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-m", "--meta"), default="NA", help="meta file information for each cell. 1 clolumn should have Cell"),
  make_option(c("--color_gene"), default="NA", help="target gene for coloring."),
  make_option(c("--sufix"), default="", help="sufix of image file name"),
  make_option(c("--color_by"), default="NA", help="color_by"),
  make_option(c("--shape_by"), default="NA", help="shape_by"),
  make_option(c("--size_by"), default="NA", help="size_by"),
  make_option(c("--color"), default="NA", help="custome colors. separated by ,"),
  make_option(c("--title"), default="NA", help="title of each graph"),
  make_option(c("--legend"), default="NA", help="legend image file"),
  make_option(c("--gene"), default="NA", help="file of gene list"),
  make_option(c("--cell"), default="NA", help="file of cell list"),
  make_option(c("--le_height"), default="NA", help="legend height"),
  make_option(c("--le_width"), default="NA", help="legend width")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mat <- as.character(opt["in"])
DIR_out <- as.character(opt["out"])
if(substring(DIR_out, nchar(DIR_out), nchar(DIR_out)) != "/"){
  DIR_out <- paste0(DIR_out, "/")
}
FILE_meta <- as.character(opt["meta"])
NAME_sufix <- as.character(opt["sufix"])
if(NAME_sufix != ""){
  NAME_sufix <- paste0("_", NAME_sufix)
}
FILE_serat <- paste0(DIR_out, "serat.rds")
TARGET_GENE <- as.character(opt["color_gene"])


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

### get coloring parameters
setOption <- function(string){
  ooo <- as.character(opt[string])
  if(ooo == "NA"){
    NULL
  }else{
    ooo
  }
}
setOption2 <- function(string){
  if(as.character(opt[string]) == "NA"){
    NULL
  }else{
    as.numeric(as.character(opt[string]))
  }
}
color_by <- setOption("color_by")
shape_by <- setOption("shape_by")
size_by <- setOption("size_by")
TITLE <- setOption("title")
FILE_legend <- setOption("legend")
le_height <- setOption2("le_height")
le_width <- setOption2("le_width")

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
}else if(TARGET_GENE != "NA"){
 if(TARGET_GENE %in% rownames(serat@scale.data)){
   D_cell <- data.frame(Cell=colnames(serat@scale.data), Val=serat@scale.data[TARGET_GENE,], stringsAsFactors = FALSE)
   colnames(D_cell)[2] <- TARGET_GENE
 }else{
   cat(TARGET, "is not found in matrices\n")
   q()
 }
}else{
  D_cell <- NULL
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
plot_PCA(pca, file=paste0(DIR_out, "pca", NAME_sufix, ".png"), cell_table = D_cell, color_by=color_by, size_by=size_by, shape_by = shape_by, pallete = Colors, title = TITLE)

if(NAME_sufix == ""){
  ok_cell <- Filter_PCA(pca, file=paste0(DIR_out, "pca_outlier.png"))
}else{
  ok_cell <- Filter_PCA(pca)
}
FILE_pca2 <- paste0(DIR_out, "pca2.rds")
if(!file.exists(FILE_pca2)){
  pca2 <- Clustering_PCA(serat, genes=Gene_list, cells=ok_cell)
  saveRDS(pca2, FILE_pca2)
}else{
  pca2 <- readRDS(FILE_pca2)
}
p1 <- plot_PCA(pca2, file=paste0(DIR_out, "pca2", NAME_sufix, ".png"), cell_table = D_cell, color_by=color_by, size_by=size_by, shape_by = shape_by, pallete = Colors, title = TITLE)



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
p2 <- plot_tSNE(tsne, file=paste0(DIR_out, "tsne", NAME_sufix, ".png"), cell_table = D_cell, color_by=color_by, size_by=size_by, shape_by = shape_by, pallete = Colors, title = TITLE, legend_file = FILE_legend, le_height = le_height, le_width = le_width)



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
p3 <- plot_SIMILR(sim, file=paste0(DIR_out, "SIMILR", NAME_sufix, ".png"), cell_table = D_cell, color_by=color_by, size_by=size_by, shape_by = shape_by, pallete = Colors, title = TITLE)


#=============================================
# Mix three graph
#=============================================
p <- plot_grid(p1 + theme(legend.position="none"), p2+ theme(legend.position="none")+ labs(title=""), p3+ theme(legend.position="none")+ labs(title=""), nrow=3, align="v")
save_plot(paste0(DIR_out, "clustering_mix", NAME_sufix, ".png"), p, ncol=1, nrow=3)



