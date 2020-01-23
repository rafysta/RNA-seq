# RNA-seq anlaysis using DESeq2


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="database file name"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-g", "--groups"), default="NA", help="group name separated by commma"),
  make_option(c("-b", "--batch"), default="NA", help="batch name separated by commma"),
  make_option(c("-s", "--samples"), default="NA", help="target column names separated by commma"),
  make_option(c("-c", "--condition"), default="all", help="condition1,condition2 for comparison (not sample name, use group name). all to output all combinations")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(cowplot)))

#=====================================
# get path of program directory
#=====================================
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)


### test data
# FILE_DB <- "T:/Project/010_20190930_Chris_RNAseq_firstTime/data/RSEM_genes_result.db"
# SAMPLES <- c("Lee_b1", "Lee_b2", "Lee_b3", "White", "WT_bio1", "WT_bio2", "Zinzen")
# group <- c("Lee", "Lee", "Lee", "White", "WT", "WT", "Zinzen")
# LIBRARY_clustering <- "T:/Library/RNA-seq/Clustering_tSNE_PCA.R"
# DIR_OUT <- "T:/Project/010_20190930_Chris_RNAseq_firstTime/out/2019-10-31_DESeq2_analysis/"


DIR_OUT <- as.character(opt["out"])
if(substring(DIR_OUT, nchar(DIR_OUT), nchar(DIR_OUT)) != "/"){
  DIR_OUT <- paste(DIR_OUT, "/", sep="")
}
SAMPLES <- unlist(strsplit(as.character(opt["samples"]), ","))
group <- unlist(strsplit(as.character(opt["groups"]), ","))



# データの読み込
FILE_DB <- as.character(opt["in"])
con = dbConnect(SQLite(), FILE_DB)
D_matrix <- dbGetQuery(con, "select * from read")
dbDisconnect(con)
rm(con)

DATA <- as.matrix(D_matrix[,-1])
rownames(DATA) <- D_matrix$gene
DATA <- DATA[,SAMPLES]

### filtering
keep <- rowSums(DATA) > 1
DATA <- DATA[keep,]

condition <- factor(group)
OPT_batch <- as.character(opt["batch"])
if(OPT_batch != "NA"){
  batch <- unlist(strsplit(OPT_batch, ","))
  batch <- factor(batch)
}
rm(D_matrix)

FILE_DESeq_result <- paste0(DIR_OUT, "deseq2.rds")
if(file.exists(FILE_DESeq_result)){
  dds <- readRDS(FILE_DESeq_result)
}else{
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  if(OPT_batch != "NA"){
    coldata <- data.frame(row.names=colnames(DATA), condition, batch)
  }else{
    coldata <- data.frame(row.names=colnames(DATA), condition)
  }
  
  # print(coldata)
  dds <- DESeqDataSetFromMatrix(countData=round(DATA), colData=coldata, design= ~ condition)
  
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  # resultsNames(dds)
  saveRDS(dds, FILE_DESeq_result)
}


# Plot dispersions
FILE_dispersion <- paste0(DIR_OUT, "qc-dispersions.png")
if(!file.exists(FILE_dispersion)){
  png(FILE_dispersion, 1000, 1000, pointsize=20)
  plotDispEsts(dds,main="Dispersion plot")
  dummy <- dev.off()
}



# Regularized log transformation for clustering/heatmaps, etc
con = dbConnect(SQLite(), FILE_DB)
df <- tryCatch(dbGetQuery(con, "select * from rlog"),
   error = function(c){
     rld <- rlogTransformation(dds)
     # head(assay(rld))
     # hist(assay(rld))

     DATA.rlog <- assay(rld)
     df <- data.frame(gene=rownames(DATA.rlog), DATA.rlog)
     dbWriteTable(con, "rlog", df, row.names= FALSE, overwrite=TRUE)
     df
   }
)
DATA.rlog <- df[,SAMPLES]
rownames(DATA.rlog) <- df$gene
df <- tryCatch(dbGetQuery(con, "select * from vsd"),
               error = function(c){
                 vsd <- vst(dds)
                 if(OPT_batch != "NA"){
                   DATA.vsd <- limma::removeBatchEffect(assay(vsd), vsd$batch)
                 }else{
                   DATA.vsd <- assay(vsd)
                 }
                 df <- data.frame(gene=rownames(DATA.vsd), DATA.vsd)
                 dbWriteTable(con, "vsd", df, row.names= FALSE, overwrite=TRUE)
                 df
               }
)
dbDisconnect(con)
DATA.vsd <- df[,SAMPLES]
rownames(DATA.vsd) <- df$gene
rm(df, con)



# Sample distance heatmap
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
FILE_distance <- paste0(DIR_OUT, "qc-heatmap-samples.png")
if(!file.exists(FILE_distance)){
  sampleDists <- as.matrix(dist(t(DATA.vsd)))
  png(FILE_distance, w=1000, h=1000, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(15, 15))
  dummy <- dev.off()
  plot(as.matrix(sampleDists))
}


if(!file.exists(paste0(DIR_OUT, "qc-pca12.png"))){
  LIBRARY_clustering <- paste(script.basename, "../../Clustering_tSNE_PCA.R", sep="/")
  FLAG_loaded <- TRUE
  source(LIBRARY_clustering)
  FILE_pca <- paste0(DIR_OUT, "pca.rds")
  if(file.exists(FILE_pca)){
    pca <- readRDS(FILE_pca)
  }else{
    pca <- Clustering_PCA(DATA.vsd)
    saveRDS(pca, FILE_pca)
  }
  width <- 4 + max(nchar(group))/8
  cell_table <- data.frame(Cell=SAMPLES, group=group, sample=SAMPLES, stringsAsFactors = FALSE)
  option <- theme(text = element_text(size=15)) + theme_bw()
  plot_PCA(pca, file=paste0(DIR_OUT, "qc-pca12.png"), cell_table=cell_table, color_by="group", size_by=1.5, label="sample", width=width, option=option)
  plot_PCA(pca, file=paste0(DIR_OUT, "qc-pca23.png"), cell_table=cell_table, color_by="group", size_by=1.5, label="sample", Xcom=2, Ycom=3, width=width, option=option)
  plot_PCA(pca, file=paste0(DIR_OUT, "qc-pca13.png"), cell_table=cell_table, color_by="group", size_by=1.5, label="sample", Xcom=1, Ycom=3, width=width, option=option)
}


# differential expression genes
getDiffGenes <- function(con1, con2){
  CON_TYPE <- paste0("DE_", con1, "_divide_", con2)
  dir.create(paste0(DIR_OUT, CON_TYPE), showWarnings = FALSE)
  DIR_out_sub <- paste0(DIR_OUT, CON_TYPE, "/")
  
  FILE_res <- paste0(DIR_out_sub, "res.rds")
  
  if(file.exists(FILE_res)){
    res <- readRDS(FILE_res)
  }else{
    res <- results(dds, contrast=c("condition", con1, con2), alpha=0.05)
    saveRDS(res, FILE_res)
  }
  # mcols(res, use.names=TRUE)
  # summary(res)
  # table(res$padj<0.05)
  
  ### Examine plot of p-values
  FILE_MA <- paste0(DIR_out_sub, "pvaldis.png")
  png(FILE_MA, 1000, 1000, pointsize=20)
  hist(res$pvalue, breaks=50, col="grey", main="Histogram of P-value", xlab="P-value")
  dummy <- dev.off()

  index_up <- res$log2FoldChange > 0 & res$padj < 0.05 & !is.na(res$padj)
  index_down <- res$log2FoldChange < 0 & res$padj < 0.05 & !is.na(res$padj)
  sig_cate <- rep(0, nrow(DATA))
  sig_cate[index_up] <- 1
  sig_cate[index_down] <- -1
  
  ### MA plot
  FILE_MA <- paste0(DIR_out_sub, "maplot.png")
  df <- data.frame(x=res$baseMean, y=res$log2FoldChange, p=res$padj)
  p <- ggplot(df, aes(x, y)) + geom_point(alpha=0.3, size=1.3, col=ifelse(!is.na(df$p) & df$p < 0.05, "red", "grey30")) +
    scale_x_log10(  breaks = trans_breaks("log10", function(x) 10^x),　labels = trans_format("log10", math_format(10^.x))  ) +
    labs(x="Average score", y=paste0("Log2 (", con1, " / ", con2, ")"), 
         subtitle=paste0("Up: ", sum(index_up), ", Down: ", sum(index_down), ", NoChange: ", length(sig_cate) - sum(index_up) - sum(index_down))) +
    geom_hline(yintercept = 0, col="blue", lty=2, lwd=1, alpha=0.7) +
    theme_bw() + theme(text=element_text(size=12))
  save_plot(FILE_MA, p, base_height = 5, base_width = 6)
  
  
  output <- data.frame(gene=rownames(res), average=res$baseMean, log2FC=res$log2FoldChange, p=res$pvalue, FDR=res$padj, SigCate=sig_cate)
  
  con = dbConnect(SQLite(), FILE_DB)
  dbWriteTable(con, CON_TYPE, output, row.names= FALSE, overwrite=TRUE)
  dbDisconnect(con)
  
  ## Write results
  FILE_DE <- paste0(DIR_out_sub, "DE.txt")
  write.table(output, file=FILE_DE, quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)
  
  number <- c(sum(index_up), length(sig_cate) - sum(index_up) - sum(index_down), sum(index_down), length(sig_cate))
  pct <- round(number / number[4] * 100, digits = 2)
  
  SUMMARY_TABLE <- cbind(number, pct)
  rownames(SUMMARY_TABLE) <- c("Up", "NoSig", "Down", "Total")
  colnames(SUMMARY_TABLE) <- c("number", "percent")
  
  write.table(SUMMARY_TABLE, file=paste0(DIR_out_sub, "Summary_", con1, "_divide_", con2, ".txt"), 
              quote = F, sep="\t", eol="\n", row.names = T, col.names = NA)
}

COMBINATIONS <- as.character(opt["condition"])
if(COMBINATIONS == "all"){
  ccc <- combn(unique(group),2)
  for(n in 1:ncol(ccc)){
    getDiffGenes(ccc[1,n], ccc[2,n])
  }
}else{
  targetConditions <- unlist(strsplit(COMBINATIONS, ","))
  getDiffGenes(targetConditions[1], targetConditions[2])
}



