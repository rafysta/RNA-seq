#!/usr/bin/Rscript
# making database from quantitate result

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="RSEM result files. separated by ,"),
  make_option(c("-n", "--names"), default="NA", help="name of each input RSEM result"),
  make_option(c("-o", "--out"), default="NA", help="output sqlite3 database file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(dplyr)))

FILE_RSEM <- unlist(strsplit(as.character(opt["in"]), ","))
NAMES <- unlist(strsplit(as.character(opt["names"]), ","))
FILE_DB <- as.character(opt["out"])

getRSEM <- function(n){
  df <- fread(FILE_RSEM[n])
  df <- df %>% select(gene=gene_id, read=expected_count, TPM, FPKM)
  df <- df %>% mutate(sample=NAMES[n])
  df
}
D_table <- do.call(rbind, pblapply(1:length(FILE_RSEM), getRSEM))

con = dbConnect(SQLite(), FILE_DB)
for(t in c("read", "TPM", "FPKM")){
  df <- D_table[,c("gene", "sample", t)] %>% tidyr::spread(key = sample, value = !!as.name(t))
  dbWriteTable(con, t, df, row.names= FALSE, overwrite=TRUE)
}
dbDisconnect(con)



