#!/usr/bin/Rscript
# CIBERSORTの結果をまとめる

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="CIBERSORT result file"),
  make_option(c("-o", "--out"), default="NA", help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))

FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])


D_CIBER <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
cell_types <- colnames(D_CIBER)
cell_types <- cell_types[!(cell_types %in% c("P-value",	"Pearson Correlation", "RMSE"))]
cell_types_out <- data.frame(Cell=rownames(D_CIBER), cell_type=cell_types[max.col(D_CIBER[,cell_types])], pval=D_CIBER[,"P-value"], stringsAsFactors = FALSE)

write.table(cell_types_out, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


cell_types_out %>% group_by(cell_type) %>% summarise(count=n()) %>% as.data.frame()
