# create R object from text file

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="Talbe file"),
  make_option(c("-o", "--out"),help="Talbe object")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_table <- as.character(opt["in"])
FILE_object <- as.character(opt["out"])


DATA_head <- read.table(FILE_table, header=TRUE, nrows = 5, check.names = FALSE, sep="\t", quote=NULL)
classes <- sapply(DATA_head, class)
DATA_all <- read.table(FILE_table, header=TRUE, check.names = FALSE, sep="\t", colClasses = classes, quote=NULL)

# Remove duplicated gene (make unique list)
DATA_all <- DATA_all[!duplicated(DATA_all[,4]),]

rownames(DATA_all) <- DATA_all[,4]


column_transcript <- 1
DATA_all <- DATA_all[, -c(1,4)]
colnames(DATA_all) <- c("chr", "strand", "symbol")

saveRDS(DATA_all, FILE_object)
