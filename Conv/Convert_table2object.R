#!/usr/bin/Rscript
# create R object from text file

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="Talbe file"),
  make_option(c("--gene"), default="NA", help="gene name conversion table file"),
  make_option(c("--cell"), default="NA", help="cell name conversion table file or database query"),
  make_option(c("--db"), default="NA", help="database file"),
  make_option(c("-o", "--out"),help="Talbe object")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(dplyr)))

FILE_table <- as.character(opt["in"])
FILE_object <- as.character(opt["out"])
DB_sample <- as.character(opt["db"])

D_head <- read.table(FILE_table, header=TRUE, nrows = 5, check.names = FALSE, sep="\t", stringsAsFactors = FALSE)
classes <- sapply(D_head, class)
D_table <- read.table(FILE_table, header=TRUE, check.names = FALSE, sep="\t", colClasses = classes, stringsAsFactors = FALSE)
rm(D_head)

### Gene table
FILE_gene <- as.character(opt["gene"])
if(FILE_gene != "NA"){
  colnames(D_table)[1] <- "transcript"
  if(file.exists(FILE_gene)){
    D_gene <- read.table(FILE_gene, header=FALSE, sep="\t", stringsAsFactors = FALSE, col.names = c("EnsembleID", "transcript"))
  }
  D_table <- dplyr::left_join(D_table, D_gene, by="transcript") %>% select(-transcript)
  
  # Change order of columns
  sample_names <- colnames(D_table)
  D_table <- D_table[,c("EnsembleID",  sample_names[sample_names != "EnsembleID"])]
}else{
  colnames(D_table)[1] <- "EnsembleID"
}


### Cell table
# 1. Column in table (only cell for output)
# 2. new name
FILE_cell <- as.character(opt["cell"])
if(FILE_cell != "NA"){
  if(file.exists(DB_sample)){
    con = dbConnect(SQLite(), DB_sample)
    D_cell <- dbGetQuery(con, FILE_cell)
    dbDisconnect(con)
  }
  
  ### Check column is exists
  if(sum(!(D_cell[,1] %in% colnames(D_table))) > 0){
    cat("Some sample doesn't exists in table\n")
    q()
  }
  D_table <- D_table[,c("EnsembleID", D_cell[,1])]
  
  rename <- function(dat, oldnames, newnames) {
    datnames <- colnames(dat)
    datnames[which(datnames %in% oldnames)] <- newnames
    colnames(dat) <- datnames
    dat
  }
  D_table <- rename(D_table, as.character(D_cell[,1]), as.character(D_cell[,2]))
}

mat <- D_table %>% select(-EnsembleID) %>% as.matrix()
rownames(mat) <- D_table %>% pull(EnsembleID)


# at least 1 read exists
mat <- mat[apply(mat, 1, sum) > 0,]


saveRDS(mat, FILE_object)



