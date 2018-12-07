#!/usr/bin/Rscript
# Automation of clustering of single cell RNA-seq 

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="file of normalized matrices"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-m", "--meta"), default="NA", help="meta file information for each cell"),
  make_option(c("--gene"), default="NA", help="file of gene list")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mat <- as.character(opt["in"])
DIR_out <- as.character(opt["out"])
if(substring(DIR_out, nchar(DIR_out), nchar(DIR_out)) != "/"){
  DIR_out <- paste0(DIR_out, "/")
}
FILE_meta <- as.character(opt["meta"])

FILE_serat <- paste0(DIR_out, "serat.rds")

Gene_list <- read.table(as.character(opt["gene"]), header=FALSE, sep="\t", stringsAsFactors = FALSE)
Gene_list <- Gene_list[,1]

#=============================================
# Serat objectを取得
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


#=============================================
# PCA
#=============================================
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))

### PCA (1回目)
serat <- RunPCA(object = serat, pc.genes = Gene_list, do.print = FALSE, pcs.compute=20)
p <- PCAPlot(object = serat, dim.1 = 1, dim.2 = 2, pt.size=1, do.return=TRUE, no.legend =TRUE, plot.title="PCA")
save_plot(paste0(DIR_out, "PCA_1.png"), p, base_height = 5, base_width = 5)


RemoveOutlier <- function(x){
  abs(x - median(x))/mad(x, constant=1)
}
pca_mat <- as.data.frame(serat@dr$pca@cell.embeddings)
mad_PC1 <- RemoveOutlier(pca_mat$PC1) < 20
mad_PC2 <- RemoveOutlier(pca_mat$PC2) < 20
cell_name <- rownames(pca_mat)

### PCAでおかしな値を示すcellを除く
ok_cell <- cell_name[mad_PC1 & mad_PC2]
serat2 <- SubsetData(object = serat, cells.use = ok_cell)

### おかしなcellの名前をプロット
D_pca <- pca_mat %>% select(PC1, PC2) %>% tidyr::gather(key = "PC", value = "score") %>% mutate(mad=c(mad_PC1, mad_PC2)) %>%
  mutate(name=if_else(mad, "", rep(cell_name,2)))
p <- ggplot(D_pca, aes(x=PC, y=score, colour=PC, fill=PC,  label =name)) + geom_violin() +
  geom_jitter(colour=if_else(D_pca %>% pull(mad), "gray20", "purple"))+
  geom_text_repel()+
  theme(legend.position="none") + labs(x="", y="PCA score")
save_plot(paste0(DIR_out, "PCA_score_distribution.png"), p, base_height = 5, base_width = 4)


### PCA(2回目)
serat2  <- RunPCA(object = serat2, pc.genes = Gene_list, do.print = FALSE, pcs.compute=20)
p <- PCAPlot(object = serat2, dim.1 = 1, dim.2 = 2, pt.size=1, do.return=TRUE, no.legend =TRUE, plot.title="PCA")
save_plot(paste0(DIR_out, "PCA_2.png"), p, base_height = 5, base_width = 5)
P1 <- p


# 追加の情報を読み取る
if(FILE_meta != "NA"){
  D_cell <- read.table(FILE_meta, header=TRUE, sep="\t", stringsAsFactors = FALSE)
}

#=============================================
# tSNE
#=============================================
serat <- RunTSNE(object = serat, dims.use = 1:10, do.fast = TRUE)
p <- TSNEPlot(object = serat, pt.size=1, do.return=TRUE, no.legend =TRUE, plot.title="tSNE")
save_plot(paste0(DIR_out, "tSNE_1.png"), p, base_height = 5, base_width = 5)

serat2 <- RunTSNE(object = serat2, dims.use = 1:10, do.fast = TRUE)
p <- TSNEPlot(object = serat2, pt.size=1, do.return=TRUE, no.legend =TRUE, plot.title="tSNE")
save_plot(paste0(DIR_out, "tSNE_2.png"), p, base_height = 5, base_width = 5)
p2 <- p

#=============================================
# SIMILR
#=============================================
suppressWarnings(suppressMessages(library(SIMLR)))
suppressWarnings(suppressMessages(library(igraph)))

pca_data <- t(serat2@dr$pca@cell.embeddings)

### 上位3位のclusterを出す
SIM_estimate <- SIMLR_Estimate_Number_of_Clusters(pca_data, 2:10, cores.ratio = 0)
SIM_estimate <- data.frame(cluster_num=2:10, value=SIM_estimate$K1)
SIM_estimate <- SIM_estimate %>% filter(min_rank(value) < 4) %>% arrange(desc(value)) %>% pull(cluster_num)

set.seed(1111)
for(cluster in SIM_estimate){
  sim <- SIMLR(X = pca_data, c = cluster)
  colnames(sim) <- c("SIMLR_1",  "SIMLR_2")
  p <- ggplot(sim, aes(x=SIMLR_1, y=SIMLR_2)) + geom_point(alpha=0.5, size=1) +
    labs(title=paste0("SIMLR (cluster# = ",cluster, ")")) + theme(legend.position="none")
  save_plot(paste0(DIR_out, "SIMLR_cluster_", cluster, ".png"), p, base_height = 5, base_width = 5)
  p3 <- p
}

p <- plot_grid(p1, p2, p3, nrow=3, align="v")
save_plot(paste0(DIR_out, "clustering_mix.png"), p, ncol=1, nrow=3)


