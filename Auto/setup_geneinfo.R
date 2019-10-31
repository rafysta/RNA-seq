#!/usr/bin/Rscript
# making database of gene information

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-x", "--organism"), default="NA", help="organism name"),
  make_option(c("-o", "--out"), default="NA", help="output sqlite3 database file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(dplyr)))

HOME <- path.expand("~")
switch (as.character(opt["organism"]),
  "fly" = {
    FILE_GO <- paste0(HOME, "/Genome/data/drosophila_melanogaster/BDGP6.22/GO.txt")
    FILE_GENE.INFO <- paste0(HOME, "/Genome/data/drosophila_melanogaster/BDGP6.22/GENE.INFO.txt")
  }
)
FILE_DB <- as.character(opt["out"])

con = dbConnect(SQLite(), FILE_DB)

GENE.INFO <- fread(FILE_GENE.INFO)
dbWriteTable(con, "gene", GENE.INFO, row.names= FALSE, overwrite=TRUE)
rm(GENE.INFO)

GO <- fread(FILE_GO)
GO.map <- GO %>% filter(accession != "") %>% select(gene, accession) %>% distinct()
dbWriteTable(con, "GO_map", GO.map, row.names= FALSE, overwrite=TRUE)
rm(GO.map)

GO.def <- GO %>% select(accession, domain, name, definition) %>% distinct()
dbWriteTable(con, "GO_def", GO.def, row.names= FALSE, overwrite=TRUE)
rm(GO.def)

dbDisconnect(con)



