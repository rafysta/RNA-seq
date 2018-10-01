# create R object from text file

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="Talbe file"),
  make_option(c("-o", "--out"),help="Talbe object")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_table <- as.character(opt["in"])
FILE_object <- as.character(opt["out"])

DATA_head <- read.table(FILE_table, header=TRUE, nrows = 5, check.names = FALSE, sep="\t")
classes <- sapply(DATA_head, class)
DATA_all <- read.table(FILE_table, header=TRUE, check.names = FALSE, sep="\t", colClasses = classes)

# 重複する遺伝子名があったら合計する
DATA_all <- aggregate(DATA_all[,-1], list(DATA_all[,1]), sum)

DATA <- DATA_all[,2:ncol(DATA_all)]
rownames(DATA) <- DATA_all[,1]


samples <- colnames(DATA)
# nameList <- strsplit(samples, ":")
# nameMatrix <- matrix(unlist(nameList), ncol=2, byrow=TRUE)


DATA <- DATA[,samples]
# at least 1 read exists
DATA <- DATA[apply(DATA, 1, sum) > 0,]


saveRDS(DATA, FILE_object)
