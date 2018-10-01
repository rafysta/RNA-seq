#!/usr/bin/Rscript
# Ontology analysis using GO slim

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), help="file with list of target genes"),
  make_option(c("-o", "--out"), help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

DB_hg19 <- "W:/Data/Human_seq/hg19/hg19.db"
DB_hg19 <- "/wistar/noma/Data/Human_seq/hg19/hg19.db"

TOTAL_gene <- 20000

FILE_in <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-06-12_GO_slim_analysis_of_compartment_switched_gene/gene_BA.txt"

FILE_in <- as.character(opt["in"])
DATA <- read.table(FILE_in, header=FALSE, as.is = TRUE, stringsAsFactors = FALSE)
genes <- DATA[,1]


library("RSQLite")
con = dbConnect(SQLite(), DB_hg19)
GOs <- dbGetQuery(con, "select distinct(acc) as acc, name, count(EnsembleID) as count from go_slim left join go_term using (acc) where hasID_in_gene=1 group by acc having count > 5;")
rownames(GOs) <- GOs[,1]

TABLE <- c()
for(acc in GOs[,1]){
  query <- paste("select EnsembleID from go_slim where acc = '", acc, "' and hasID_in_gene = 1;", sep="")
  mother <- dbGetQuery(con, query)
  mother <- mother[,1]
  
  a <- sum(genes %in% mother)
  b <- length(genes) - a
  c <- length(mother) - a
  d <- TOTAL_gene - a - b - c
  
  mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
  ff <- fisher.test(mx, alternative = "greater")
  Pval <- format(ff$p.value, digits = 3)
  Enrich <- as.numeric(format(ff$estimate, digits = 3))
  result <- c(GOs[acc,"acc"], GOs[acc,"name"], GOs[acc,"count"], a, Enrich, Pval)
  TABLE <- rbind(TABLE, result)
}
dbDisconnect(con)
colnames(TABLE) <- c("Accession numbers", "GO names", "Total genes", paste("Gene in target(n=", length(genes), ")", sep=""), "Enrichment", "P-values")


write.table(TABLE, as.character(opt["out"]), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


