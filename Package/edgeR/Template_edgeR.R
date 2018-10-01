# EdgeRによる解析


#======================================================
# 以下の環境を設定する
#======================================================

# プロジェクト名
PROJECT_NAME <- "2017-02-23_RNA-seq_human_IMR90_BJ_Growing_Senescent"

# 入出力フォルダ名
INPUT_OUTPUT_DIRECTORY <- "2017-02-27_significantly_different_expressed_genes"

# sampleのグループ設定
group <- c("G1", "G1", "G2", "S1", "S1", "S2", "G3", "G3", "G3", "S3", "S3", "S3")

# もしbiological replicaが存在しない場合、以下の２つカラムのデータをreplicaとみなす
fake_replica <- c(2, 3)

# ここまで
#======================================================


### データの読み込み
FILE_rawdata <- paste("W:/Project/", PROJECT_NAME, "/out/", INPUT_OUTPUT_DIRECTORY, "/count_all.txt", sep="")
DIR_out <-paste("W:/Project/", PROJECT_NAME, "/out/", INPUT_OUTPUT_DIRECTORY, "/", sep="")

rawdata <- read.table(FILE_rawdata, check.names=FALSE, stringsAsFactors = FALSE, sep="\t", header=TRUE)

# sampleの数
sampleNum <- ncol(rawdata)-1
sampleName <- colnames(rawdata)[-1]

if(sampleNum != length(group)){
  cat ("Please check group setting")
  q()
}



library(edgeR)
y <- DGEList(counts=rawdata[,2:(sampleNum+1)], genes=rawdata[,1], group=group)


# design matrixの作成
design <- model.matrix(~group)


# lowly expressed geneを除く
# keep <- rowSums(cpm(y) > 1) >= 2
keep <- rowSums(cpm(y) > 1) >= 10
y <- y[keep, keep.lib.sizes=FALSE]


# TMM normalizationをする
y <- calcNormFactors(y)
y$samples

plotMDS(y)

#======================================================
# Dispersionの計算
#======================================================
if(length(group) == length(unique(group))){
  # Biological replicaが無い場合の計算
  # Dispersionを推測する（ばらつき具合） XXXとXXXのデータについては特に変化が無いと仮定
  y1 <- DGEList(counts=rawdata[,fake_replica], genes=rawdata[,1], group=c(1,1))
  y0 <- estimateCommonDisp(y1)
  y$common.dispersion <- y0$common.dispersion
}else{
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
}


# 異なって発現する遺伝子を見つける
DifferentExpressedGene <- function(group1, group2){
  et <- exactTest(y, pair=c(group1, group2)) # group1とgroup2 を比較する
  
  
  # すべての結果をファイルに出力する
  FILE_out <- paste(DIR_out, "comparison_table_", group1, "_vs_", group2, ".txt", sep="")
  out <- topTags(et, n = nrow(y$counts))
  write.table(file=FILE_out, x = out,sep = "\t", eol = "\n",quote = FALSE, col.names = TRUE, row.names = FALSE)


  # FDR < 0.05のものをsignificantとする
  de <- decideTestsDGE(et, adjust.method="BH", p.value=0.05)
  index_significant <- rownames(y)[as.logical(de)]
    
  # average logCPMとlogFCのプロット
  FILE_out <- paste(DIR_out, "graph_", group1, "_vs_", group2, ".png", sep="")
  png(filename=FILE_out, width=600, height=500)
  par(oma=c(0,0,0,0), mar=c(4,4,3,1))
  plotSmear(et, main=paste("logFC ( ", group2, " / ", group1, " )", sep=""), cex=0.35, de.tags=index_significant)
  abline(h = c(-2, 2), col = "blue")
  dev.off()
    
  summary(de)
}

# group 
group

DifferentExpressedGene("G1", "S1")
DifferentExpressedGene("G2", "S2")
DifferentExpressedGene("G3", "S3")








