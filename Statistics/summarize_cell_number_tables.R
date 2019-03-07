#!/usr/bin/Rscript
# cell numberについての表を３つ出力する

#=============================================
# 設定
#=============================================
DIR <- "X:/hideki_projects/423_20190226_singleCell_CTCL/"
FILE_DB <- paste0(DIR, "db/Data.db")
FILE_rsem <- paste0(DIR, "data/rsem_count.txt")

DIR_out <- paste0(DIR, "out/2019-03-04_make_several_graphs_for_quality_check/")

#=============================================

suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggsci)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(ggrepel)))


### Control以外のcellのbarcode
con = dbConnect(SQLite(), FILE_DB)
D_table <- dbGetQuery(con, "select * from well")
dbDisconnect(con)
rm(con)
barcode_not_control <- D_table %>% filter(!grepl("Ctrl", Sample2)) %>% pull(Cell)
rm(D_table)


### RSEMをmatrixに変換
D_rsem <- read.table(FILE_rsem, header=TRUE, sep="\t", stringsAsFactors = FALSE)
D_rsem <- D_rsem %>% filter(Cell %in% barcode_not_control)
DATA <- D_rsem %>% select(Gene, Cell, Count) %>% tidyr::spread(key = Cell, value = Count, fill=0)
rm(D_rsem)
rownames(DATA) <- DATA[,1]
DATA <- DATA[,-1]
DATA <- as.matrix(DATA)


lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}


# at least M% cells have N genes detected with >=T counts
# cat("Calc at least M% cells have N genes detected with >=T counts...", format(Sys.time(), "%X"), "\n")
X_scales <- c(1, 2, 5, 10, 20, 50, 100, 500, 1000, 5000, 10000)
Num_more_T <- c()
for(T in X_scales){
  scores <- c()
  Passed_GeneNum_per_cell <- sort(apply(DATA >= T, 2, sum), decreasing = TRUE)
  for(M in c(5, 25, 50, 75, 95)){
    scores <- c(scores, as.numeric(Passed_GeneNum_per_cell[length(Passed_GeneNum_per_cell)*M/100]))
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- paste(c(5, 25, 50, 75, 95), "%", sep="")
rownames(Num_more_T) <- X_scales
write.table(Num_more_T, paste0(DIR_out, "Gene_number_having_at_least_Xpercent_cell_has_Kreads.txt"), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")



# at least M% cells have common N genes detected with >=T counts
# cat("Calc at least M% cells have common N genes detected with >=T counts...", format(Sys.time(), "%X"), "\n")
Num_more_T <- c()
for(T in X_scales){
  scores <- c()
  Passed_CellNum_per_gene <- sort(apply(DATA > T, 1, sum), decreasing = TRUE)
  for(M in c(5, 25, 50, 75, 95)){
    scores <- c(scores, sum(Passed_CellNum_per_gene > ncol(DATA)*M/100))
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- paste(c(5, 25, 50, 75, 95), "%", sep="")
rownames(Num_more_T) <- X_scales
write.table(Num_more_T, paste(DIR_out, "Common_Gene_number_for_Xpercent_cell_having_K_reads.txt", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")



### cell number having N genes with more than T reads
# cat("Calc cell number having N genes with more than T reads...", format(Sys.time(), "%X"), "\n")
X_scales <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
Num_more_T <- c()
for(T in c(1, 5, 10, 50, 100)){
  scores <- c()
  Passed_geneNum_per_cell <- apply(DATA >= T, 2, sum)
  for(N in X_scales){
    scores <- c(scores, sum(Passed_geneNum_per_cell >= N))
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- X_scales
rownames(Num_more_T) <- c(1, 5, 10, 50, 100)
write.table(t(Num_more_T), paste(DIR_out, "Cell_number_having_K_genes_with_X_reads.txt", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")



