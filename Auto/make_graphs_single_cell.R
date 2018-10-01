# create R object from text file

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-t", "--table"),help="Talbe object"),
  make_option(c("-q", "--quality"),help="quality file"),
  make_option(c("--align"),help="align stat file"),
  make_option(c("--annotation"),help="annotation object file"),
  make_option(c("-o", "--out"),help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_object <- "X:/hideki_projects/364_20180215_SingleCellRNA_Human/out/2018-02-15_automation/Read_count.rds"
DIR_OUT <- "X:/hideki_projects/364_20180215_SingleCellRNA_Human/out/2018-02-15_automation/img/"
FILE_quality <- "X:/hideki_projects/364_20180215_SingleCellRNA_Human/out/2018-02-15_automation/Quality_per_cell.txt"
FILE_align <- "X:/core_project/MohsenMohamed/20180212_SingleCell/Results/AlignStat.txt"
FILE_annotation <- "X:/hideki_projects/364_20180215_SingleCellRNA_Human/out/2018-02-15_automation/Read_annotation.rds"

FILE_object <- as.character(opt["table"])
DIR_OUT <- as.character(opt["out"])
FILE_quality <- as.character(opt["quality"])
FILE_align <- as.character(opt["align"])
FILE_annotation <- as.character(opt["annotation"])
RequireRead <- 100000

colorRB <- colorRampPalette(c(rgb(85,72,193, max=255), rgb(88,76,196, max=255), rgb(90,79,199, max=255), rgb(92,83,202, max=255), rgb(94,87,205, max=255), rgb(96,90,208, max=255), rgb(98,94,211, max=255), rgb(100,97,214, max=255), rgb(103,101,216, max=255), rgb(105,104,219, max=255), rgb(107,108,221, max=255), rgb(109,111,224, max=255), rgb(111,114,226, max=255), rgb(114,118,229, max=255), rgb(116,121,231, max=255), rgb(118,124,233, max=255), rgb(120,128,235, max=255), rgb(122,131,237, max=255), rgb(125,134,239, max=255), rgb(127,137,240, max=255), rgb(129,140,242, max=255), rgb(131,144,244, max=255), rgb(134,147,245, max=255), rgb(136,150,246, max=255), rgb(138,153,248, max=255), rgb(140,156,249, max=255), rgb(143,158,250, max=255), rgb(145,161,251, max=255), rgb(147,164,252, max=255), rgb(149,167,253, max=255), rgb(152,169,253, max=255), rgb(154,172,254, max=255), rgb(156,175,254, max=255), rgb(159,177,255, max=255), rgb(161,180,255, max=255), rgb(163,182,255, max=255), rgb(165,184,255, max=255), rgb(168,187,255, max=255), rgb(170,189,255, max=255), rgb(172,191,255, max=255), rgb(174,193,255, max=255), rgb(176,195,254, max=255), rgb(179,197,254, max=255), rgb(181,199,253, max=255), rgb(183,201,253, max=255), rgb(185,203,252, max=255), rgb(187,204,251, max=255), rgb(189,206,250, max=255), rgb(192,207,249, max=255), rgb(194,209,248, max=255), rgb(196,210,246, max=255), rgb(198,211,245, max=255), rgb(200,213,244, max=255), rgb(202,214,242, max=255), rgb(204,215,240, max=255), rgb(206,216,239, max=255), rgb(208,217,237, max=255), rgb(209,218,235, max=255), rgb(211,218,233, max=255), rgb(213,219,231, max=255), rgb(215,219,229, max=255), rgb(217,220,227, max=255), rgb(218,220,224, max=255), rgb(220,221,222, max=255), rgb(222,220,219, max=255), rgb(224,219,216, max=255), rgb(225,218,214, max=255), rgb(227,217,211, max=255), rgb(229,216,208, max=255), rgb(230,215,205, max=255), rgb(232,213,202, max=255), rgb(233,212,199, max=255), rgb(234,210,196, max=255), rgb(236,209,193, max=255), rgb(237,207,190, max=255), rgb(238,205,187, max=255), rgb(239,203,183, max=255), rgb(239,201,180, max=255), rgb(240,199,177, max=255), rgb(241,197,174, max=255), rgb(242,195,171, max=255), rgb(242,193,168, max=255), rgb(242,191,165, max=255), rgb(243,188,161, max=255), rgb(243,186,158, max=255), rgb(243,183,155, max=255), rgb(243,181,152, max=255), rgb(243,178,149, max=255), rgb(243,175,146, max=255), rgb(243,173,142, max=255), rgb(243,170,139, max=255), rgb(242,167,136, max=255), rgb(242,164,133, max=255), rgb(241,161,130, max=255), rgb(241,158,127, max=255), rgb(240,155,124, max=255), rgb(239,152,121, max=255), rgb(239,148,118, max=255), rgb(238,145,115, max=255), rgb(237,142,111, max=255), rgb(235,138,109, max=255), rgb(234,135,106, max=255), rgb(233,131,103, max=255), rgb(232,128,100, max=255), rgb(230,124,97, max=255), rgb(229,120,94, max=255), rgb(227,117,91, max=255), rgb(225,113,88, max=255), rgb(224,109,85, max=255), rgb(222,105,83, max=255), rgb(220,101,80, max=255), rgb(218,97,77, max=255), rgb(216,93,75, max=255), rgb(214,89,72, max=255), rgb(212,84,69, max=255), rgb(209,80,67, max=255), rgb(207,76,64, max=255), rgb(205,71,62, max=255), rgb(202,66,59, max=255), rgb(200,61,57, max=255), rgb(197,56,55, max=255), rgb(194,51,52, max=255), rgb(192,45,50, max=255), rgb(189,39,48, max=255), rgb(186,33,46, max=255), 
                              rgb(183,25,44, max=255), rgb(180,15,41, max=255), rgb(177,1,39, max=255)), space="Lab")( 100 )


ANNO <- readRDS(FILE_annotation)
DATA <- readRDS(FILE_object)
nameList <- strsplit(colnames(DATA), ":")
nameMatrix <- matrix(unlist(nameList), ncol=2, byrow=TRUE)
rownames(nameMatrix) <- colnames(DATA)
table(as.character(nameMatrix[,1]))



Quality <- as.matrix(read.table(FILE_quality, header=TRUE, check.names=FALSE, sep = "\t", row.names=2))
BarcodeList <- rownames(Quality)
Category_quality <- as.character(unique(Quality[,"State"]))
colorCategory_quality <- c("#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", 
                     "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080")[1:length(Category_quality)]
colorCells_quality <- rep('gray', nrow(Quality))
for(i in 1:length(Category_quality)){
  index <- which(as.character(Quality[,"State"]) == Category_quality[i])
  colorCells_quality[index] <- colorCategory_quality[i]
}

Category_sample <- as.character(unique(Quality[,"sample"]))
colorCategory_sample <- c("#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", 
                   "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080")[1:length(Category_sample)]
colorCells_sample <- rep('gray', nrow(Quality))
for(i in 1:length(Category_sample)){
  index <- which(as.character(Quality[,"sample"]) == Category_sample[i])
  colorCells_sample[index] <- colorCategory_sample[i]
}


# Normalize read number for each cells.
TotalRead_per_cells <- apply(DATA, 2, sum)
names(TotalRead_per_cells) <- as.character(nameMatrix[,2])
DATA.normalized <- t(t(DATA) / TotalRead_per_cells * 1e6)
DATA.normalized.log2 <- log2(DATA.normalized + 1)


# read数の色
LogTotalRead_per_cells <- log(TotalRead_per_cells)
min <- min(LogTotalRead_per_cells); max <- max(LogTotalRead_per_cells)
colReadAbs <- as.integer((LogTotalRead_per_cells - min)/(max-min)*100)
colReadAbs[which(colReadAbs==0)] = 1
colReadAbs <- colorRB[colReadAbs]


# sample name
FILE_OUT <- paste(DIR_OUT, "category_quality_name_vertical.png", sep="")
png(file=FILE_OUT, width=max(nchar(Category_quality))*10 + 10, height=18*length(Category_quality), units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(rep(2, length(Category_quality)),  length(Category_quality):1, col=colorCategory_quality, pch=15, cex=2,
     axes=F, xlab="", ylab="", xlim=c(1, max(nchar(Category_quality))*7 + 3), bty = "n", ylim=c(0.5, length(Category_quality)+0.5))
text(rep(7, length(Category_quality)), length(Category_quality):1, labels=Category_quality, cex= 1, pos=4)
dummy<-dev.off()

FILE_OUT <- paste(DIR_OUT, "category_quality_name_horizontal.png", sep="")
png(file=FILE_OUT, width=18*length(Category_quality), height=max(nchar(Category_quality))*7 + 10, units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(1:length(Category_quality), rep(2, length(Category_quality)),   col=colorCategory_quality, pch=20, cex=2,
     axes=F, xlab="", ylab="", xlim=c(0.5, length(Category_quality)+0.5), bty = "n", ylim=c(1, max(nchar(Category_quality))*7 + 3))
text(1:length(Category_quality), rep(10, length(Category_quality)), labels=Category_quality, cex= 1, pos=4, srt=90, offset = -0.02)
dummy<-dev.off()

FILE_OUT <- paste(DIR_OUT, "category_sample_name_vertical.png", sep="")
png(file=FILE_OUT, width=max(nchar(Category_sample))*10 + 10, height=18*length(Category_sample), units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(rep(2, length(Category_sample)),  length(Category_sample):1, col=colorCategory_sample, pch=15, cex=2,
     axes=F, xlab="", ylab="", xlim=c(1, max(nchar(Category_sample))*7 + 3), bty = "n", ylim=c(0.5, length(Category_sample)+0.5))
text(rep(7, length(Category_sample)), length(Category_sample):1, labels=Category_sample, cex= 1, pos=4)
dummy<-dev.off()

FILE_OUT <- paste(DIR_OUT, "category_sample_name_horizontal.png", sep="")
png(file=FILE_OUT, width=18*length(Category_sample), height=max(nchar(Category_sample))*7 + 10, units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(1:length(Category_sample), rep(2, length(Category_sample)),   col=colorCategory_sample, pch=20, cex=2,
     axes=F, xlab="", ylab="", xlim=c(0.5, length(Category_sample)+0.5), bty = "n", ylim=c(1, max(nchar(Category_sample))*7 + 3))
text(1:length(Category_sample), rep(10, length(Category_sample)), labels=Category_sample, cex= 1, pos=4, srt=90, offset = -0.02)
dummy<-dev.off()


###################################################
# read distribution for each cells
###################################################
png(filename=paste(DIR_OUT, "read_number_for_each_cells_quality.png", sep=""), width=12.6, height=12.6, units = "cm", res=72)
par(mar=c(1,5,1,1), oma=c(0,0,0,0))
plot(TotalRead_per_cells[BarcodeList], log='y', pch=20, cex=1, xlab="", ylab="Total read",col=colorCells_quality,cex.lab = 1.5,xaxt="n")
abline(h=RequireRead, col='blue', lty=2, lwd=2)
mtext(as.character(RequireRead), side=2, outer=FALSE, at=RequireRead, line = 1, font=2)
dummy<-dev.off()

png(filename=paste(DIR_OUT, "read_number_for_each_cells_sample.png", sep=""), width=12.6, height=12.6, units = "cm", res=72)
par(mar=c(1,5,1,1), oma=c(0,0,0,0))
plot(TotalRead_per_cells[BarcodeList], log='y', pch=20, cex=1, xlab="", ylab="Total read",col=colorCells_sample,cex.lab = 1.5,xaxt="n")
abline(h=RequireRead, col='blue', lty=2, lwd=2)
mtext(as.character(RequireRead), side=2, outer=FALSE, at=RequireRead, line = 1, font=2)
dummy<-dev.off()


Cell_more_than_requiredRead_quality <- c()
for(S in Category_quality){
  Cell_more_than_requiredRead_quality <- c(Cell_more_than_requiredRead_quality, sum(TotalRead_per_cells[rownames(Quality[which(Quality[,"State"]==S),])] > RequireRead))
}
names(Cell_more_than_requiredRead_quality) <- Category_quality
T <- table(Quality[,"State"])
DF <- cbind(Cell_more_than_requiredRead_quality, T[Category_quality])
colnames(DF) <- c(paste(">", RequireRead, sep=""), "Total")
write.table(DF, paste(DIR_OUT, "read_number_more_than_threshold_quality.txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE, col.names = NA)


Cell_more_than_requiredRead_sample <- c()
for(S in Category_sample){
  Cell_more_than_requiredRead_sample <- c(Cell_more_than_requiredRead_sample, sum(TotalRead_per_cells[rownames(Quality[which(Quality[,"sample"]==S),])] > RequireRead))
}
names(Cell_more_than_requiredRead_sample) <- Category_sample
T <- table(Quality[,"sample"])
DF <- cbind(Cell_more_than_requiredRead_sample, T[Category_sample])
colnames(DF) <- c(paste(">", RequireRead, sep=""), "Total")
write.table(DF, paste(DIR_OUT, "read_number_more_than_threshold_sample.txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE, col.names = NA)

index_control <- as.character(nameMatrix[,1]) == "Pos" | as.character(nameMatrix[,1]) == "Neg"
DF <- data.frame(category=as.character(nameMatrix[index_control,1]), barcode=as.character(nameMatrix[index_control,2]), total_read=TotalRead_per_cells[index_control])
write.table(DF, paste(DIR_OUT, "read_number_more_than_threshold_for_control.txt", sep=""), quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE)



###################################################
# accumulated read number per cell
###################################################
cat("Calc accumulated read number per cell...", format(Sys.time(), "%X"), "\n")
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}
x_scales <- lseq(1, max(TotalRead_per_cells), length.out = 100)
Num_more_T <- list()
for(T in x_scales){
  for(S in Category_quality){
    Num_more_T[[S]] <- c(Num_more_T[[S]], sum(TotalRead_per_cells[rownames(Quality[which(Quality[,"State"]==S),])] > T))
  }
}
ylim=c(0,max(unlist(Num_more_T)))
FILE_png <- paste(DIR_OUT, "read_number_for_each_cell_accumulated_threshold_quality.png", sep="")
png(filename=FILE_png, width=12.6, height=9, units = "cm", res=72)
par(mar=c(4,4,2,1), oma=c(0,0,0,0))
plot(x_scales, Num_more_T[[Category_quality[1]]], log='x', type='l', ylab='Cell number', cex.lab = 1.4, col=colorCategory_quality[1], xlab='Read threshold', lwd=2, ylim=ylim)
for(i in 2:length(Category_quality)){
  par(new=T)
  plot(x_scales, Num_more_T[[Category_quality[i]]], log='x', type='l', axes=F, ylab='', xlab='', col=colorCategory_quality[i], lwd=2, ylim=ylim, cex.lab = 1.4)
}
abline(v=RequireRead, col='blue', lty=2, lwd=2, untf=TRUE)
mtext(as.character(RequireRead), side=3, outer=FALSE, at=RequireRead, line = 0.3, font=2)
dummy <- dev.off()

Num_more_T <- list()
for(T in x_scales){
  for(S in Category_sample){
    Num_more_T[[S]] <- c(Num_more_T[[S]], sum(TotalRead_per_cells[rownames(Quality[which(Quality[,"sample"]==S),])] > T))
  }
}
ylim=c(0,max(unlist(Num_more_T)))
FILE_png <- paste(DIR_OUT, "read_number_for_each_cell_accumulated_threshold_sample.png", sep="")
png(filename=FILE_png, width=12.6, height=9, units = "cm", res=72)
par(mar=c(4,4,2,1), oma=c(0,0,0,0))
plot(x_scales, Num_more_T[[Category_sample[1]]], log='x', type='l', ylab='Cell number', cex.lab = 1.4, col=colorCategory_sample[1], xlab='Read threshold', lwd=2, ylim=ylim)
for(i in 2:length(Category_sample)){
  par(new=T)
  plot(x_scales, Num_more_T[[Category_sample[i]]], log='x', type='l', axes=F, ylab='', xlab='', col=colorCategory_sample[i], lwd=2, ylim=ylim, cex.lab = 1.4)
}
abline(v=RequireRead, col='blue', lty=2, lwd=2, untf=TRUE)
mtext(as.character(RequireRead), side=3, outer=FALSE, at=RequireRead, line = 0.3, font=2)
dummy <- dev.off()



###################################################
# gene # more than threshold
###################################################
cat("Calc gene # more than threshold...", format(Sys.time(), "%X"), "\n")
MaxRead_per_gene <- apply(DATA[,!index_control], 1, max)

lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}
X_scales <- c(lseq(1, 10000, length.out = 10))


Num_more_T <- c()
for(T in X_scales){
  Num_more_T <- c(Num_more_T, sum(MaxRead_per_gene >= T))
}

ymax <- max(Num_more_T)
FILE_png <- paste(DIR_OUT, "Read_number_vs_gene_number.png", sep="")
png(filename=FILE_png, width=20, height=12.6, units = "cm", res=72)
par(mar=c(4,4.5,1,1), oma=c(0,0,0,0))
plot(X_scales, Num_more_T, type='o', ylab='# of gene', cex=2.2, cex.lab = 1.8, cex.axis = 1.5, xlab='Read# threshold', lwd=2, log='xy', ylim=c(1,ymax), pch=20)
text(X_scales, Num_more_T, labels=Num_more_T, cex= 1.5, pos=1)
dummy <- dev.off()


# at least M% cells have N genes detected with >=T counts
cat("Calc at least M% cells have N genes detected with >=T counts...", format(Sys.time(), "%X"), "\n")
X_scales <- c(1, 2, 5, 10, 20, 50, 100, 500, 1000, 5000, 10000)
Num_more_T <- c()
for(T in X_scales){
  scores <- c()
  Passed_GeneNum_per_cell <- sort(apply(DATA[,!index_control] >= T, 2, sum), decreasing = TRUE)
  for(M in c(5, 25, 50, 75, 95)){
    scores <- c(scores, Passed_GeneNum_per_cell[length(Passed_GeneNum_per_cell)*M/100])
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- paste(c(5, 25, 50, 75, 95), "%", sep="")
rownames(Num_more_T) <- X_scales
write.table(Num_more_T, paste(DIR_OUT, "Gene_number_having_at_least_Xpercent_cell_has_Kreads.txt", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")

# at least M% cells have common N genes detected with >=T counts
cat("Calc at least M% cells have common N genes detected with >=T counts...", format(Sys.time(), "%X"), "\n")
Num_more_T <- c()
for(T in X_scales){
  scores <- c()
  Passed_CellNum_per_gene <- sort(apply(DATA[,!index_control] > T, 1, sum), decreasing = TRUE)
  for(M in c(5, 25, 50, 75, 95)){
    scores <- c(scores, sum(Passed_CellNum_per_gene > ncol(DATA[,!index_control])*M/100))
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- paste(c(5, 25, 50, 75, 95), "%", sep="")
rownames(Num_more_T) <- X_scales
write.table(Num_more_T, paste(DIR_OUT, "Common_Gene_number_for_Xpercent_cell_having_K_reads.txt", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")



### cell number having N genes with more than T reads
cat("Calc cell number having N genes with more than T reads...", format(Sys.time(), "%X"), "\n")
X_scales <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
Num_more_T <- c()
for(T in c(1, 5, 10, 50, 100)){
  scores <- c()
  Passed_geneNum_per_cell <- apply(DATA[,nameMatrix[,2] %in% BarcodeList] >= T, 2, sum)
  for(N in X_scales){
    scores <- c(scores, sum(Passed_geneNum_per_cell >= N))
  }
  Num_more_T <- rbind(Num_more_T, scores)
}
colnames(Num_more_T) <- X_scales
rownames(Num_more_T) <- c(1, 5, 10, 50, 100)
write.table(t(Num_more_T), paste(DIR_OUT, "Cell_number_having_K_genes_with_X_reads.txt", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")





###################################################
# sample quality
###################################################
cat("Calc sample quality...", format(Sys.time(), "%X"), "\n")

for(q in c("Signal", "IntegSignal","Circularity", "Size", "Confidence", "DropIndex")){
  ymin <- min(as.numeric(Quality[BarcodeList,q]), na.rm=TRUE)
  ymax <- max(as.numeric(Quality[BarcodeList,q]), na.rm=TRUE)
  
  yrange <- c(ymin, ymax) # 範囲
  tickValues <- pretty(yrange) # 自動で区分け
  if(ymin > 1000){
    tickStrings <- format(tickValues, scientific = TRUE)
  }else{
    tickStrings <- tickValues
  }
  
  
  FILE_png <- paste(DIR_OUT, "Quality_", q, "_quality.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  plot(Quality[BarcodeList,q], col=colorCells_quality, pch=20, cex=1, cex.main = 1, xlab="", xaxt="n", ylab="",
       ylim=yrange, yaxt="n", main=q)
  axis(side=2,at=tickValues,labels=tickStrings, las=2, cex.axis=1)
  dummy <- dev.off()
  
  FILE_png <- paste(DIR_OUT, "Quality_", q, "_sample.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  plot(Quality[BarcodeList,q], col=colorCells_sample, pch=20, cex=1, cex.main = 1, xlab="", xaxt="n", ylab="",
       ylim=yrange, yaxt="n", main=q)
  axis(side=2,at=tickValues,labels=tickStrings, las=2, cex.axis=1)
  dummy <- dev.off()
  
  
  DATA_for_boxplot <- list()
  for(m in Category_sample){
    index <- Quality[,"sample"] == m
    DATA_for_boxplot[[m]] <- as.numeric(Quality[BarcodeList[index],q])
    DATA_for_boxplot[[m]] <- ifelse(is.na(DATA_for_boxplot[[m]]), ymin, DATA_for_boxplot[[m]])
  }
  FILE_png <- paste(DIR_OUT, "Quality_box_", q, "_sample.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  boxplot(DATA_for_boxplot, las = 2, col=colorCategory_sample, ylab="", main=q, cex.main = 1, ylim=yrange, yaxt="n")
  axis(side=2,at=tickValues,labels=tickStrings, las=2, cex.axis=1)
  dummy <- dev.off()
  
  
  DATA_for_boxplot <- list()
  for(m in Category_quality){
    index <- Quality[,"State"] == m
    DATA_for_boxplot[[m]] <- as.numeric(Quality[BarcodeList[index],q])
    DATA_for_boxplot[[m]] <- ifelse(is.na(DATA_for_boxplot[[m]]), ymin, DATA_for_boxplot[[m]])
  }
  FILE_png <- paste(DIR_OUT, "Quality_box_", q, "_quality.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  boxplot(DATA_for_boxplot, las = 2, col=colorCategory_quality, ylab="", main=q, cex.main = 1, ylim=yrange, yaxt="n")
  axis(side=2,at=tickValues,labels=tickStrings, las=2, cex.axis=1)
  dummy <- dev.off()
}




###################################################
# alignment status
###################################################
cat("Calc alignment status...", format(Sys.time(), "%X"), "\n")
AlignStat <- as.matrix(read.table(FILE_align, header=TRUE, check.names=FALSE, row.names=1))
for(q in c("Contam", "Transcript", "Genome", "RSEM")){
  
  q1 <- paste("Unique", q, sep="")
  q2 <- paste("Multi", q, sep="")
  
  if(q == "Contam"){
    q = "Ribosomal RNA"
  }
  
  
  scores <- (as.numeric(AlignStat[BarcodeList,q1]) + as.numeric(AlignStat[BarcodeList,q2]))/ as.numeric(AlignStat[BarcodeList,"Total"]) * 100
  
  FILE_png <- paste(DIR_OUT, "Alignment_", q, "_quality.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  plot(scores, col=colorCells_quality, pch=20, cex=1, cex.lab = 1.5, xlab="", xaxt="n", ylim=c(0,100), ylab="", 
       main=paste(q, "(%)"), cex.main = 1, las=2, cex.axis=1)
  dummy <- dev.off()
  
  FILE_png <- paste(DIR_OUT, "Alignment_", q, "_sample.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  plot(scores, col=colorCells_sample, pch=20, cex=1, cex.lab = 1.5, xlab="", xaxt="n", ylim=c(0,100), ylab="", 
       main=paste(q, "(%)"), cex.main = 1, las=2, cex.axis=1)
  dummy <- dev.off()
  
  
  DATA_for_boxplot <- list()
  for(m in Category_quality){
    index <- Quality[,"State"] == m
    DATA_for_boxplot[[m]] <- (as.numeric(AlignStat[BarcodeList[index],q1]) + as.numeric(AlignStat[BarcodeList[index],q2]))/ as.numeric(AlignStat[BarcodeList[index],"Total"]) * 100
  }
  FILE_png <- paste(DIR_OUT, "Alignment_box_", q, "_quality.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  boxplot(DATA_for_boxplot, las = 2, col=colorCategory_quality, ylab="", main=paste(q, "(%)"), cex.main = 1, ylim=c(0,100), cex.axis=1)
  dummy <- dev.off()
  
  DATA_for_boxplot <- list()
  for(m in Category_sample){
    index <- Quality[,"sample"] == m
    DATA_for_boxplot[[m]] <- (as.numeric(AlignStat[BarcodeList[index],q1]) + as.numeric(AlignStat[BarcodeList[index],q2]))/ as.numeric(AlignStat[BarcodeList[index],"Total"]) * 100
  }
  FILE_png <- paste(DIR_OUT, "Alignment_box_", q, "_sample.png", sep="")
  png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
  par(mar=c(4,4,2,1), oma=c(0,0,0,0))
  boxplot(DATA_for_boxplot, las = 2, col=colorCategory_sample, ylab="", main=paste(q, "(%)"), cex.main = 1, ylim=c(0,100), cex.axis=1)
  dummy <- dev.off()
}




###################################################
# % of mitochondoria 
# GRCh38にMitochondoriaが入っていないので飛ばす
###################################################
# cat("Calc % a of mitochondoria ...", format(Sys.time(), "%X"), "\n")
# {
#   index_Mitochondoria_gene <- ANNO[rownames(DATA), "chr"] == "MT"
#   TotalMitochondorial_read_per_cells <- apply(DATA[index_Mitochondoria_gene,], 2, sum)
#   names(TotalMitochondorial_read_per_cells) <- nameMatrix[,2]
#   scores <- TotalMitochondorial_read_per_cells / TotalRead_per_cells * 100
#   names(scores) <- as.character(nameMatrix[,2])
#   
#   FILE_png <- paste(DIR_OUT, "Alignment_mitochondoria.png", sep="")
#   png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
#   par(mar=c(4,4,2,1), oma=c(0,0,0,0))
#   plot(scores[BarcodeList], col=colorCells, pch=20, cex=1, cex.lab = 1.5, xlab="", xaxt="n", ylim=c(0,100), ylab="", 
#        main="Mitochondoria Transcript(%)", cex.main = 1, las=2, cex.axis=1)
#   dummy <- dev.off()
#   
#   
#   OUTPUT_DATA <- cbind(CategoryTable[BarcodeList,1], TotalRead_per_cells[BarcodeList], TotalMitochondorial_read_per_cells[BarcodeList])
#   colnames(OUTPUT_DATA) <- c("Category", "Total_transcript", "Mitochondoria_aligned_transcript")
#   write.table(OUTPUT_DATA, paste(DIR_OUT, "Summary_of_mitochondoria_aligned_reads.txt", sep=""), 
#               col.names=NA, row.names = TRUE, sep="\t", quote = FALSE, eol = "\n")
#   rm(OUTPUT_DATA)
#   
#   DATA_for_boxplot <- list()
#   for(m in Category){
#     index <- CategoryTable == m
#     DATA_for_boxplot[[m]] <- scores[BarcodeList[index]]
#   }
#   
#   FILE_png <- paste(DIR_OUT, "Alignment_box_mitochondoria.png", sep="")
#   png(filename=FILE_png, width=8, height=8, units = "cm", res=72)
#   par(mar=c(4,4,2,1), oma=c(0,0,0,0))
#   boxplot(DATA_for_boxplot, las = 2, col=colorCategory, ylab="", main="Mitochondoria Transcript(%)", cex.main = 1, ylim=c(0,100), cex.axis=1)
#   dummy <- dev.off()
# }


###################################################
# actinの数
###################################################
index_actin <- rownames(ANNO[ANNO[,"symbol"] == "ACTB", ])
actin_threshold <- qnorm(0.01, mean=mean(DATA.normalized.log2[index_actin,]), sd=sd(DATA.normalized.log2[index_actin,]))

# 各cellのactinの数
FILE_png <- paste(DIR_OUT, "Actin_read_distribution.png", sep="")
png(filename=FILE_png, width=12, height=12, units = "cm", res=72)
par(mar=c(4,4,2,1), oma=c(0,0,0,0))
plot(DATA.normalized.log2[index_actin,], pch=20, cex=1, cex.lab = 1.5, xlab="", ylab="Normalized score of Actin", cex.axis=1, col=colReadAbs,xaxt="n")
abline(h=actin_threshold, col='green', lty=2)
text(1, actin_threshold, labels = "0.01% line", pos = 4)
dummy <- dev.off()

# Scaleの出力
FILE_png <- paste(DIR_OUT, "Color_scale.png", sep="")
png(filename=FILE_png, width=4, height=1.5, units = "cm", res=72, bg="transparent")
color_data <- matrix(1:length(colorRB), ncol=1)
haba <- length(colorRB)/4
par(oma=c(0,0,0,0), mar=c(1,0,0,0))
image(color_data, col=colorRB, axes=F)
label.loc <- c(1, haba, 2*haba, 3*haba, length(colorRB))/100
axis(1, line = 0.2, at=label.loc, cex.axis=1.4)
dummy <- dev.off()


# nomralized read < 2 のtotal cell の分布
FILE_png <- paste(DIR_OUT, "Actin_low_read_cell_hist.png", sep="")
png(filename=FILE_png, width=12, height=10, units = "cm", res=72, bg="transparent")
par(oma=c(0,0,0,0), mar=c(4,4,2,1))
hist(TotalRead_per_cells, col=rgb(0.9,0,0,0.5),  breaks = seq(0,max(TotalRead_per_cells), length.out = 100), 
     xlim=c(0,4e6), xlab="Total read# / cell", main="", cex.lab=1.4, cex.axis=1.2)
hist(TotalRead_per_cells[nameMatrix[colnames(DATA.normalized.log2[,DATA.normalized.log2[index_actin,] < 2]),2]], 
     col=rgb(0,0,0.9,0.5), add=T, breaks = seq(0,max(TotalRead_per_cells), length.out = 100), xlim=c(0,4e6))
abline(v=RequireRead, col='green', lty=2, lwd=2.5)
mtext(as.character(RequireRead), side=3, outer=FALSE, at=RequireRead, line = 0.3, font=2, cex=1.4)
dummy <- dev.off()


sum(TotalRead_per_cells[nameMatrix[colnames(DATA.normalized.log2[,DATA.normalized.log2[index_actin,] < 2]),2]] > RequireRead)



###################################################
# Gene saturationの計算
###################################################
cat("Gene saturation calculation ...", format(Sys.time(), "%X"), "\n")
{
  ### 3個のcellで9以上のreadがある遺伝子のみをターゲットにする
  DATA2 <- DATA[apply(DATA[,!index_control]>9, 1, sum) > 3,]
  colnames(DATA2) <- nameMatrix[,2]
  
  randSampling <- function(ori, num){
    ori <- ori[ori > 0]
    csum <- cumsum(ori)
    TotalRead <- sum(ori)
    if(num > TotalRead){
      NA
    }else if(num == TotalRead){
      ori
    }else if(num > TotalRead/2){
      cate <- cut(sample(TotalRead, TotalRead - num, replace = FALSE), c(0, csum))
      result <- ori - table(cate)
      as.numeric(result)
    }else{
      cate <- cut(sample(TotalRead, num, replace = FALSE), c(0, csum))
      as.numeric(table(cate))
    }
  }
  
  readNumbers <- seq(4e4, 1e6, by=4e4)
  Barcode_GOOD_and_notControl <- intersect(rownames(Quality[Quality[,"State"] == "Good",]), nameMatrix[!index_control,2])
  Max_read_barcode <- names(TotalRead_per_cells[TotalRead_per_cells==max(TotalRead_per_cells[Barcode_GOOD_and_notControl])])
  
  SATURATION <- c()
  for(Threshold_min_read in c(0, 1, 2, 5, 10, 20,50,100, 200)){
    GeneNumbers <- c()
    for(T in readNumbers){
      sampled <- randSampling(DATA2[,Max_read_barcode], T)
      GeneNumbers <- c(GeneNumbers, sum(sampled > Threshold_min_read))
    }
    SATURATION <- cbind(SATURATION, GeneNumbers)
  }
  
  Threshold <- c(0, 1, 2, 5, 10, 20,50,100, 200)
  colnames(SATURATION) <- Threshold
  rownames(SATURATION) <- readNumbers
  write.table(SATURATION, file=paste(DIR_OUT, "saturation_table_raw.txt", sep=""), row.names = TRUE, col.names = NA, sep="\t", eol = "\n", quote = FALSE)
  
  col <- rainbow(length(Threshold))
  for(i in 1:length(Threshold)){
    ylim <- c(0, max(SATURATION[,i]) * 1.2)
    png(paste(DIR_OUT, "saturation_plot_read_", Threshold[i], ".png", sep=""), width=12.6, height=13, units = "cm", res=72)
    par(oma=c(0,0,0,0), mar=c(4,5,3,1))
    plot(as.numeric(rownames(SATURATION))/1000, SATURATION[,i], type='b', col=col[i], lwd=1.5, ylim=ylim,
         xlab="Reads (x1000)", pch=20, cex=1.2, cex.lab=2, cex.axis=1.8, ylab=paste("Gene # (reads > ", Threshold[i], ")", sep=""))
    abline(v=RequireRead/1000, col='blue', lty=2, lwd=2)
    mtext(as.character(RequireRead/1000), side=3, outer=FALSE, at=RequireRead/1000, line = 0.3, font=2, cex=2)
    dev.off()
  }
}


###################################################
# Z scoreの計算
###################################################
cat("Calc correlations ...", format(Sys.time(), "%X"), "\n")
if(file.exists(paste(DIR_OUT, "Correlation_matrix.rds", sep=""))){
  cor.matrix <- readRDS(paste(DIR_OUT, "Correlation_matrix.rds", sep=""))
}else{
  calcCor <- function(i){
    cor(DATA[,i], DATA, method = "spearman")
  }
  cor.matrix <- cbind(sapply(1:ncol(DATA), calcCor))
  colnames(cor.matrix) <- colnames(DATA)
  rownames(cor.matrix) <- colnames(DATA)
  saveRDS(cor.matrix, paste(DIR_OUT, "Correlation_matrix.rds", sep=""))
}

# Zscore calculation
cat("Calc Z score ...", format(Sys.time(), "%X"), "\n")
{
  if(file.exists(paste(DIR_OUT, "Zscore.rds", sep=""))){
    Zscore <- readRDS(paste(DIR_OUT, "Zscore.rds", sep=""))
  }else{
    calcZscore <- function(i){
      s_med <- median(as.numeric(cor.matrix[i,]))
      t_med <- median(as.numeric(cor.matrix))
      sd_mad <- mad(as.numeric(cor.matrix))
      (s_med - t_med)/sd_mad
    }
    Zscore <- sapply(1:nrow(cor.matrix), calcZscore)
    names(Zscore) <- colnames(DATA)
    saveRDS(Zscore, paste(DIR_OUT, "Zscore.rds", sep=""))
  }
  # read数とZ-scoreの関係
  FILE_png <- paste(DIR_OUT, "Zscore_vs_readNum_quality.png", sep="")
  png(filename=FILE_png, width=15, height=10, units = "cm", res=72)
  par(mar=c(4,5,2,1), oma=c(0,0,0,0))
  plot(TotalRead_per_cells, Zscore, pch=20, cex=1, cex.lab = 1.5, xlab="Total read / cell", ylab="Z-score", 
       cex.main = 1, cex.axis=1, col=colorCells_quality, log="x")
  abline(v=RequireRead, col='blue', lty=2, lwd=2)
  mtext(as.character(RequireRead), side=3, outer=FALSE, at=RequireRead, line = 0.3, font=2)
  dummy <- dev.off()
  
  FILE_png <- paste(DIR_OUT, "Zscore_vs_readNum_sample.png", sep="")
  png(filename=FILE_png, width=15, height=10, units = "cm", res=72)
  par(mar=c(4,5,2,1), oma=c(0,0,0,0))
  plot(TotalRead_per_cells, Zscore, pch=20, cex=1, cex.lab = 1.5, xlab="Total read / cell", ylab="Z-score", 
       cex.main = 1, cex.axis=1, col=colorCells_sample, log="x")
  abline(v=RequireRead, col='blue', lty=2, lwd=2)
  mtext(as.character(RequireRead), side=3, outer=FALSE, at=RequireRead, line = 0.3, font=2)
  dummy <- dev.off()
  
  # 各cellのZ-score
  FILE_png <- paste(DIR_OUT, "Zscore_for_each_cell.png", sep="")
  png(filename=FILE_png, width=12, height=10, units = "cm", res=72)
  par(mar=c(1,5,1,1), oma=c(0,0,0,0))
  plot(Zscore, pch=20, cex=1, cex.lab = 1.5, xlab="", ylab="Z-score", cex.axis=1, col=colReadAbs,xaxt="n")
  dummy <- dev.off()
}




