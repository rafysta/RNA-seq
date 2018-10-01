# RNA-seq anlaysis using DESeq2

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="object or raw read matrices"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-g", "--groups"), default="NA", help="group name separated by commma"),
  make_option(c("-s", "--samples"), default="NA", help="target column names separated by commma"),
  make_option(c("-c", "--condition"), default="NA", help="condition1,condition2 for comparison (not sample name, use group name)"),
  make_option(c("-p", "--figure"), default=FALSE, help="output figure (TRUE) or not")
)
opt <- parse_args(OptionParser(option_list=option_list))

#=====================================
# get path of program directory
#=====================================
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

FLAG_cairo <- FALSE
if(suppressPackageStartupMessages(require(Cairo))){
  suppressPackageStartupMessages(library(Cairo))
  FLAG_cairo <- TRUE
}


DIR_OUT <- as.character(opt["out"])
if(substring(DIR_OUT, nchar(DIR_OUT), nchar(DIR_OUT)) != "/"){
  DIR_OUT <- paste(DIR_OUT, "/", sep="")
}
SAMPLES <- unlist(strsplit(as.character(opt["samples"]), ","))
group <- unlist(strsplit(as.character(opt["groups"]), ","))
FLAG_figure <- eval(parse(text=as.character(opt["figure"])))

# データの読み込み
FILE_mat <- as.character(opt["in"])
DATA <- readRDS(FILE_mat)
DATA <- DATA[,SAMPLES]
condition <- factor(group)

FILE_DESeq_result <- paste(DIR_OUT, "deseq2.rds", sep="")
if(file.exists(FILE_DESeq_result)){
  dds <- readRDS(FILE_DESeq_result)
}else{
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  coldata <- data.frame(row.names=colnames(DATA), condition)
  dds <- DESeqDataSetFromMatrix(countData=round(DATA), colData=coldata, design=~condition)
  
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  saveRDS(dds, FILE_DESeq_result)
}


# Plot dispersions
FILE_dispersion <- paste(DIR_OUT, "qc-dispersions.png", sep="")
if(FLAG_figure){
  if(FLAG_cairo){
    Cairo(file=FILE_dispersion, 1000, 1000, pointsize=20)
  }else{
    png(FILE_dispersion, 1000, 1000, pointsize=20)
  }
  plotDispEsts(dds,main="Dispersion plot")
  dummy <- dev.off()
}



# Regularized log transformation for clustering/heatmaps, etc
FILE_table_normalized_object <- paste(DIR_OUT, "read_table_norm.rds", sep="")
FILE_talbe_normalized_text <- paste(DIR_OUT, "read_table_norm.txt" ,sep="")
if(file.exists(FILE_table_normalized_object)){
  DATA.norm <- readRDS(FILE_table_normalized_object)
}else{
  rld <- rlogTransformation(dds)
  # head(assay(rld))
  # hist(assay(rld))
  DATA.norm <- assay(rld)
  saveRDS(DATA.norm, FILE_table_normalized_object)
  write.table(DATA.norm, FILE_talbe_normalized_text, sep="\t", quote = F, col.names=NA, row.name=T)
}

library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]

# Sample distance heatmap

FILE_distance <- paste(DIR_OUT, "qc-heatmap-samples.png", sep="")

if(FLAG_figure){
  suppressPackageStartupMessages(library(gplots))
  rld <- rlogTransformation(dds)
  sampleDists <- as.matrix(dist(t(DATA.norm)))
  if(FLAG_cairo){
    Cairo(file=FILE_distance, 1000, 1000, pointsize=20)
  }else{
    png(FILE_distance, w=1000, h=1000, pointsize=20)
  }
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(15, 15), main="Sample Distance Matrix")
  dummy <- dev.off()
}

if(FLAG_figure){
  LIBRARY_clustering <- paste(script.basename, "Clustering_tSNE_PCA.R", sep="/")
  FLAG_loaded <- TRUE
  source(LIBRARY_clustering)
  colors_pca <- rainbow(length(SAMPLES))
  pca <- Clustering_PCA(DATA.norm, rownames(DATA.norm), colnames(DATA.norm))
  plot_PCA(pca, col_list=colors_pca, col_eachcell=colors_pca, names=SAMPLES, file=paste(DIR_OUT, "qc-pca12.png", sep=""))
  plot_PCA(pca, col_list=colors_pca, col_eachcell=colors_pca, names=SAMPLES, file=paste(DIR_OUT, "qc-pca23.png", sep=""), Xcom=2, Ycom=3)
  plot_PCA(pca, col_list=colors_pca, col_eachcell=colors_pca, names=SAMPLES, file=paste(DIR_OUT, "qc-pca13.png", sep=""), Xcom=1, Ycom=3)
}

# differential expression genes
getDiffGenes <- function(con1, con2){
  dir.create(paste(DIR_OUT, "DE_", con1, "_", con2, sep=""), showWarnings = FALSE)
  DIR_out_sub <- paste(DIR_OUT, "DE_", con1, "_", con2, "/", sep="")
  
  FILE_res <- paste(DIR_out_sub, "res.rds", sep="")
  
  if(file.exists(FILE_res)){
    res <- readRDS(FILE_res)
  }else{
    res <- results(dds, contrast=c("condition", con1, con2), alpha=0.05)
    saveRDS(res, FILE_res)
  }
  mcols(res, use.names=TRUE)
  # summary(res)
  # table(res$padj<0.05)
  
  ### Examine plot of p-values
  FILE_MA <- paste(DIR_out_sub, "pvaldis.png", sep="")
  if(FLAG_figure){
    if(FLAG_cairo){
      Cairo(file=FILE_MA, 1000, 1000, pointsize=20)
    }else{
      png(FILE_MA, 1000, 1000, pointsize=20)
    }
    hist(res$pvalue, breaks=50, col="grey")
    dummy <- dev.off()
  }
  
  ### MA plot
  FILE_MA <- paste(DIR_out_sub, "maplot.png", sep="")
  if(FLAG_figure){
    if(FLAG_cairo){
      Cairo(file=FILE_MA, width=15, height=13, units = "cm", res = 72)
    }else{
      png(FILE_MA, width=15, height=13, units = "cm", res = 72)
    }
    ymax <- max(res$log2FoldChange, na.rm = TRUE)
    ymin <- min(res$log2FoldChange, na.rm = TRUE)
    plotMA(res, ylim=c(ymin, ymax), alpha=0.05, main="", xlab="Average", ylab=paste("log2 (", con1, "/", con2, ")", sep=""), cex.lab=1.4, cex.axis=1.2)
    # plotMA(res, alpha=0.05, main="", xlab="Average", ylab=paste("log2 (", con1, "/", con2, ")", sep=""), cex.lab=1.4, cex.axis=1.2)
    dummy <- dev.off()
  }
  
  
  index_up <- res$log2FoldChange > 0 & res$padj < 0.05 & !is.na(res$padj)
  index_down <- res$log2FoldChange < 0 & res$padj < 0.05 & !is.na(res$padj)
  sig_cate <- rep(0, nrow(DATA))
  sig_cate[index_up] <- 1
  sig_cate[index_down] <- -1
  
  output <- cbind(rownames(DATA), res$baseMean, res$log2FoldChange, res$pvalue, res$padj, sig_cate)
  colnames(output) <- c("ID", "average", "log2FC", "p", "FDR", "sigCate")
  
  ## Write results
  FILE_DE <- paste(DIR_out_sub, "DE.txt", sep="")
  write.table(output, file=FILE_DE, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)
  
  number <- c(sum(index_up), length(sig_cate) - sum(index_up) - sum(index_down), sum(index_down), length(sig_cate))
  pct <- round(number / number[4] * 100, digits = 2)
  
  SUMMARY_TABLE <- cbind(number, pct)
  rownames(SUMMARY_TABLE) <- c("Up", "NoSig", "Down", "Total")
  colnames(SUMMARY_TABLE) <- c("number", "percent")
  
  write.table(SUMMARY_TABLE, file=paste(DIR_out_sub, "Summary_", con1, "_divide_", con2, ".txt", sep=""), 
              quote = F, sep="\t", eol="\n", row.names = T, col.names = NA)
}



targetConditions <- unlist(strsplit(as.character(opt["condition"]), ","))

getDiffGenes(targetConditions[1], targetConditions[2])







  
