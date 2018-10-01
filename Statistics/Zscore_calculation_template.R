# Z-score を計算するためのプログラムのテンプレート


DIR_OUT <- "X:/hideki_projects/328_20170810_SingleCellRNA_Human/out/2017-08-24_Zscore/"

###################################################
# DATAの準備
###################################################
FILE_object <- "X:/hideki_projects/328_20170810_SingleCellRNA_Human/out/2017-08-10_automation/Read_count.rds"

DATA <- readRDS(FILE_object)
nameList <- strsplit(colnames(DATA), ":")
nameMatrix <- matrix(unlist(nameList), ncol=2, byrow=TRUE)
colnames(DATA) <- nameMatrix[,2]

# Normalize read number for each cells.
Barcode_GOOD_cell <- names(Cell_status[Cell_status == "good"])
TotalRead_per_cells <- apply(DATA, 2, sum)
AverageRead_per_cells <- mean(TotalRead_per_cells[nameMatrix[,2] %in% Barcode_GOOD_cell])
Scale_factor <- AverageRead_per_cells / TotalRead_per_cells
DATA.normalized <- t(t(DATA) * Scale_factor)
DATA.normalized.log2 <- log2(DATA.normalized + 1)

# Data for analysis
DATA_for_analysis <- DATA.normalized.log2[, Barcode_GOOD_cell]


  


###################################################
# Z scoreの計算
###################################################
# Correlation calculations
calcCor <- function(i){
  cor(DATA_for_analysis[,i], DATA_for_analysis, method = "spearman")
}
if(file.exists(paste(DIR_OUT, "Correlation_matrix.rds", sep=""))){
  cor.matrix <- readRDS(paste(DIR_OUT, "Correlation_matrix.rds", sep=""))
}else{
  cor.matrix <- cbind(sapply(1:ncol(DATA_for_analysis), calcCor))
  saveRDS(cor.matrix, paste(DIR_OUT, "Correlation_matrix.rds", sep=""))
}

# Zscore calculation
calcZscore <- function(i){
  s_med <- median(as.numeric(cor.matrix[i,]))
  t_med <- median(as.numeric(cor.matrix))
  sd_mad <- mad(as.numeric(cor.matrix))
  (s_med - t_med)/sd_mad
}

Zscore <- sapply(1:nrow(cor.matrix), calcZscore)
saveRDS(Zscore, paste(DIR_OUT, "Zscore.rds", sep=""))









