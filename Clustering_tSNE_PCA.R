#!/usr/bin/Rscript
# tSNE and PCA clustering library

if(!exists("FLAG_loaded")){
  suppressPackageStartupMessages(library("optparse"))
  option_list <- list(
    make_option(c("-r", "--read"), default="NA", help="normalized read"),
    make_option(c("-d", "--db"), default="NA", help="SQLite3 database"),
    make_option(c("-i", "--in"), default="NA", help="if tsne or pca was performed before, specify the rds file"),
    make_option(c("-o", "--out"), default="NA", help="output file prefix"),
    make_option(c("--gene"), default="NA", help="SQLite3 query to select id list or file name of gene list(first column is gene)"),
    make_option(c("--cell"), default="NA", help="SQLite3 query to output Barcode list or file name of gene list(first column is cell)"),
    make_option(c("--cell_table"), default="NA", help="SQLite3 query to output cell_table or Robject .rds file of cell_table dataframe"),
    make_option(c("-m", "--method"), default="NA", help="tSNE or PCA. Method to apply "),
    make_option(c("--perplexity"), default=10, help="In case of tSNE, which perlexity will be check?"),
    make_option(c("--graph"), default=FALSE, help="plot graph (TRUE) or not (FALSE). Default: FALSE"),
    make_option(c("--color_by"), default="NA", help="output pic were color by"),
    make_option(c("--shape_by"), default="NA", help="output pic were shape by"),
    make_option(c("--size_by"), default="NA", help="output pic were size by"),
    make_option(c("--pallete"), default="NA", help="specify colors"),
    make_option(c("--width"), default=5.2, help="width of output image")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
}


checkParameter <- function(name){
  if(as.character(opt[name]) == "NA"){
    cat("Parameter", name, "is not specified\n")
    q()
  }
}


Clustering_PCA <- function(DATA){
  DATA <- as.matrix(DATA)
  storage.mode(DATA) <- "double"
  
  var_line <- apply(DATA, 1, stats::var, na.rm = TRUE)
  keep <- is.finite(var_line) & (var_line > 0)
  if (any(!keep)) DATA <- DATA[keep, , drop = FALSE]
  
  stats::prcomp(t(DATA), center = TRUE, scale. = TRUE)
}

Clustering_tSNE <- function(DATA, seed=20, perplexity = 10, theta = 0, max_iter=3000){
  library(Rtsne)
  set.seed(seed) # 再現性の確保
  tsne <- Rtsne(t(DATA), dims = 2, verbose = FALSE, perplexity = perplexity, 
                theta = theta, max_iter=max_iter)
  tsne[["cell"]] <- colnames(DATA)
  tsne
}

plot_tSNE <- function(tsne, file="", title=NULL,  cell_table="", color_by=NULL, shape_by=NULL, size_by=NULL, 
                      width=5.5, alpha=0.5, pallete = NULL, option=NULL){
  ### cell_tableはdata.frame
  # cellというカラムが定義されていること！
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  D_score <- data.frame(Cell=tsne$cell, x=tsne$Y[,1], y=tsne$Y[,2], stringsAsFactors = FALSE)
  D_table <- dplyr::left_join(D_score, cell_table, by="Cell", copy=FALSE)
  if(is.null(color_by)){
    p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
      labs(x="Component1", y="Component2", title=title)
  }else{
    if(length(unique(D_table[,color_by])) < 20){
      if(is.null(pallete)){
        pallete <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      }
      colors <- colorRampPalette(pallete)(length(unique(D_table[,color_by])))
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=factor(D_table[,color_by]), size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_manual(values=colors, name=color_by) +
        labs(x="Component1", y="Component2", title=title)
    }else{
      mid <- median(D_table[,color_by], na.rm = TRUE)
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0)) + geom_point(alpha=alpha) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="grey90", high="red", space ="Lab" )+
        labs(x="Component1", y="Component2", title=title)
    }
  }
  if(!is.null(option)){
    p <- p + option
  }
  if(file != ""){
    save_plot(file, p, base_height = 4, base_width = width)
  }else{
    print(p)
  }
}

plot_PCA <- function(pca, file="", title=NULL, Xcom=1, Ycom=2, cell_table="", color_by=NULL, shape_by=NULL, 
                     size_by=NULL, width=5.3, alpha=0.5, label = NULL, pallete = NULL,  option=NULL){
  ### cell_tableはdata.frame
  # cellというカラムが定義されていること！
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  suppressWarnings(suppressMessages(library(ggrepel)))
  D_score <- data.frame(Cell=rownames(pca$x), x=pca$x[,Xcom], y=pca$x[,Ycom], stringsAsFactors = FALSE)
  D_table <- dplyr::left_join(D_score, cell_table, by="Cell", copy=FALSE)

  # 寄与率
  contribution <- pca$sdev^2/sum(pca$sdev^2)*100
  
  if(is.null(color_by)){
    p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0, label=label)) + geom_point(alpha=alpha) +
      labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
           y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
  }else{
    if(length(unique(D_table[,color_by])) < 20){
      if(is.null(pallete)){
        pallete <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
      }
      colors <- colorRampPalette(pallete)(length(unique(D_table[,color_by])))
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=factor(D_table[,color_by]), size=size_by, shape=shape_by, stroke=0, label=label)) + geom_point(alpha=alpha) +
        scale_color_manual(values=colors, name=color_by) + guides(size = "none") +
        labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
             y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
    }else{
      mid <- median(D_table[,color_by], na.rm = TRUE)
      p <- ggplot(D_table, aes_string(x="x", y="y", colour=color_by, size=size_by, shape=shape_by, stroke=0, label=label)) + geom_point(alpha=alpha) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="grey90", high="red", space ="Lab" )+
        labs(x=paste("PC", Xcom, " (", format(contribution[Xcom], digits = 3), "%)", sep=""),
             y=paste("PC", Ycom, " (", format(contribution[Ycom], digits = 3), "%)", sep=""), title=title)
    }
  }
  if(!is.null(label)){
    p <- p + geom_text_repel(col='black', size=2)
  }
  if(!is.null(option)){
    p <- p + option
  }
  if(file != ""){
    save_plot(file, p, base_height = 4, base_width = width)
  }else{
    print(p)
  }
}



if(!exists("FLAG_loaded")){
  library(RSQLite)
  checkParameter("read")
  checkParameter("out")
  checkParameter("method")
  
  FILE_normalized <- as.character(opt["read"])
  FILE_DB <- as.character(opt["db"])
  DIR_OUT <- as.character(opt["out"])
  METHOD <- as.character(opt["method"])
  FLAG_plot <- eval(parse(text=as.character(opt["graph"])))
  img_width <- as.numeric(as.character(opt["width"]))
  
  pallete <- as.character(opt["pallete"])
  if(pallete != "NA"){
    pallete <- strsplit(pallete, ",")
  }
  
  NA_or_character <- function(name){
    ttt = as.character(opt[name])
    if(ttt == "NA"){
      NULL
    }else{
      ttt
    }
  }
  color_by=NA_or_character("color_by")
  shape_by=NA_or_character("shape_by")
  size_by=NA_or_character("size_by")
  if(FLAG_plot){
    checkParameter("cell_table")
  }
  
  
  # normalized read
  DATA <- readRDS(FILE_normalized)
  
  
  ### Cell list
  if(as.character(opt["cell"]) != "NA"){
    if(file.exists(as.character(opt["cell"]))){
      D <- read.table(as.character(opt["cell"]), header=FALSE, sep="\t", stringsAsFactors = FALSE)
      CELLs <- D[,1]
    }else{
      library("RSQLite")
      DB <- dbConnect(SQLite(), FILE_DB)
      CELLs <- dbGetQuery(DB, as.character(opt["cell"]))[,1]
      dbDisconnect(DB)
    }
    CELLs <- intersect(CELLs, colnames(DATA))
    DATA <- DATA[,CELLs]
  }
  ### Gene list
  if(as.character(opt["gene"]) != "NA"){
    if(file.exists(as.character(opt["gene"]))){
      D <- read.table(as.character(opt["gene"]), header=FALSE, sep="\t", stringsAsFactors = FALSE)
      GENEs <- D[,1]
    }else{
      library("RSQLite")
      DB <- dbConnect(SQLite(), FILE_DB)
      GENEs <- dbGetQuery(DB, as.character(opt["gene"]))[,1]
      dbDisconnect(DB)
    }
    GENEs <- intersect(GENEs, rownames(DATA))
    DATA <- DATA[GENEs,]
  }
  ### Cell table
  if(FLAG_plot){
    if(file.exists(as.character(opt["cell_table"]))){
      cell_table <- readRDS(as.character(opt["cell_table"]))
    }else{
      library("RSQLite")
      DB <- dbConnect(SQLite(), FILE_DB)
      cell_table <- dbGetQuery(DB, as.character(opt["cell_table"]))
      dbDisconnect(DB)
    }
  }
  
  if(METHOD=="tSNE"){
    perplexity <- as.numeric(as.character(opt["perplexity"]))
    FILE_tsne <- as.character(opt["in"])
    if(FILE_tsne == "NA"){
      FILE_tsne <- paste(DIR_OUT, "tSNE_perplexity_", perplexity, ".rds", sep="")
    }
    if(file.exists(FILE_tsne)){
      tsne <- readRDS(FILE_tsne)
    }else{
      # theta=0 is original tSNE
      tsne <- Clustering_tSNE(DATA, theta = 0, max_iter = 3000, perplexity = perplexity)
      saveRDS(tsne, FILE_tsne)
    }
    if(FLAG_plot){
      FILE_png <- paste(DIR_OUT, "tSNE_perplexity_", perplexity, ".png", sep="")
      plot_tSNE(tsne, file=FILE_png, title=paste("tSNE (perplexity=", perplexity, ")", sep=""), 
                cell_table=cell_table, color_by=color_by, shape_by=shape_by, size_by=size_by, width=img_width)
    }
  }else if(METHOD == "PCA"){
    FILE_pca <- as.character(opt["in"])
    if(FILE_pca == "NA"){
      FILE_pca <- paste(DIR_OUT, "pca.rds", sep="")
    }
    if(file.exists(FILE_pca)){
      pca <- readRDS(FILE_pca)
    }else{
      pca <- Clustering_PCA(DATA)
      saveRDS(pca, FILE_pca)
    }
    if(FLAG_plot){
      xs=c(1,1,2)
      ys=c(2,3,3)
      for(i in 1:3){
        x <- xs[i]
        y <- ys[i]
        FILE_png <- paste(DIR_OUT, "PCA_PC", x, "_PC", y, ".png", sep="")
        plot_PCA(pca, file=FILE_png, title="PCA", Xcom=x, Ycom=y, cell_table=cell_table, 
                 color_by=color_by, shape_by=shape_by, size_by=size_by, width=img_width) 
      }
    }
  }
}
