 library("RColorBrewer")
 library("DESeq2")
 library("geneplotter")
 library("EDASeq")
 library("genefilter")
 library("pheatmap")
 
 setwd("/home/mirahan/Work/TE-pre/L1HSstats.DEGES")
 
 
 batchnormal<- read.table("BatchData.normal.tsv", header = TRUE, sep="\t")
 barcode <- strsplit(as.character(batchnormal[,1]), "-")
 barcode <- matrix(unlist(barcode), ncol = 7, byrow = TRUE)
 p2batchnormal <- cbind.data.frame(barcode[,3], batchnormal[,2])
 colnames(p2batchnormal) <- c("patient", "batch_normal")
 
 batchcancer<- read.table("BatchData.cancer.tsv", header = TRUE, sep="\t")
 barcode <- strsplit(as.character(batchcancer[,1]), "-")
 barcode <- matrix(unlist(barcode), ncol = 7, byrow = TRUE)
 p2batchcancer <- cbind.data.frame(barcode[,3], batchcancer[,2])
 colnames(p2batchcancer) <- c("patient", "batch_cancer")
 
 p2batch <- merge(p2batchnormal, p2batchcancer, by="patient")
 
 
 patient2cancer <- read.table("patientToCancerMappings.tsv")
 barcode <- strsplit(as.character(patient2cancer[,1]), "-")
 barcode <- matrix(unlist(barcode), ncol = 3, byrow = TRUE)
 p2ct <- cbind.data.frame(barcode[,3], patient2cancer[,2])
 colnames(p2ct) <- c("patient", "type")
 
 p2ct2batch <- merge(p2ct, p2batch, by="patient")
 #p2ct2batch <- merge(p2ct, p2batchnormal, by="patient")
 
 cancertypes = c("BLCA", "BRCA", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")
 
 if (file.exists("ddsGenes.DEGES.RData")) {
  load("ddsGenes.DEGES.RData")
 } else {
  
  file.names <- dir("RawCnts", pattern =".genes.cntTable")
  patientIDs <- substr(file.names, 9, 12)
  patient_files <- cbind.data.frame(patientIDs, file.names)
  colnames(patient_files) <- c("patient", "filename")
  patient_types <- merge(patient_files, p2ct2batch, by = "patient")
  patient_types <- patient_types[patient_types$type!="CESC", ]
  patient_types <- patient_types[patient_types$type!="PAAD", ]
  patient_types <- patient_types[patient_types$type!="PCPG", ]
  patient_types <- patient_types[patient_types$type!="SARC", ]
  patient_types <- patient_types[patient_types$type!="THYM", ]
  patient_types$type <- droplevels(patient_types$type)
  
  genenames <- rownames(read.table(paste("RawCnts",as.character(patient_types[1,2]),sep="/"), header=TRUE, row.names=1))
  #samplenames <- paste(rep(patient_types[,1], each=2), rep(c("cancer", "normal"), nrow(patient_types)), sep="_")
  #coldata <- cbind.data.frame(rep(patient_types[,1], each=2), rep(patient_types[,3], each=2), rep(c("cancer", "normal"), nrow(patient_types)), as.vector(t(patient_types[,4:5])))
  coldata <- cbind.data.frame(rep(patient_types[,1], each=1), rep(patient_types[,3], each=1), rep(c("normal"), nrow(patient_types)), as.vector(t(patient_types[,4])))
  #rownames(coldata) <- samplenames 
  rownames(coldata) <- patient_types[,1] 
  colnames(coldata) <- c("patient", "type", "condition", "batch")
  coldata$batch <- factor(coldata$batch)
  
  #need to make batch nested within type 
  coldata <- cbind.data.frame(coldata, rep(0, nrow(coldata)))
  colnames(coldata)[5] <- "nested.batch"
  for (i in 1:length(cancertypes)) {
    coldataLUAD <- coldata[coldata$type==cancertypes[i],]
    coldataLUAD$batch <- droplevels(coldataLUAD$batch)
    coldata[coldata$type==cancertypes[i],"nested.batch"]=as.numeric(coldataLUAD$batch)
  }
  coldata$nested.batch = factor(coldata$nested.batch)
  dim(coldata)
 
  cntmatrix = c() 
  if (file.exists("cntmatrix.gene.txt")) {
    cntmatrix = read.table(file="cntmatrix.gene.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
    cntGenes <- cntmatrix
  } else { 
    cntmatrix = matrix(rep(0, nrow(patient_types)*length(genenames)),  nrow=length(genenames))
    for (i in 1:nrow(patient_types)) {
      cnts <- read.table(paste("RawCnts",as.character(patient_types[i,2]),sep="/"), header=TRUE, row.names=1)
      #cntmatrix[,(i*2-1):(i*2)]=as.matrix(cnts)
      cntmatrix[,i]=as.matrix(cnts)
    }
    dim(cntmatrix)
    rownames(cntmatrix) <- genenames
    cntGenes <- cntmatrix
    colnames(cntmatrix) <- coldata$patient
    write.table(cntmatrix, file="cntmatrix.gene.txt", quote=FALSE, row.names=TRUE, sep="\t")
  } 
  
  ddsnormal <- DESeqDataSetFromMatrix(
    countData = cntGenes,
    colData = coldata,
    design = ~ nested.batch + type )
  #ddsnormal <- ddsnormal[, ddsnormal$type!="CESC"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="PAAD"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="PCPG"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="SARC"]
  #ddsnormal <- ddsnormal[, ddsnormal$type!="THYM"]
  
  ddsnormal$type <- droplevels(ddsnormal$type)
  #keep <- rowSums(counts(ddsnormal)) >= 10
  #dds <- dds[keep,]
 
  # step 1
  FDR <- 0.1
  floorPDEG <- 0.05
  filter <- apply(counts(ddsnormal), 1, function(x) { all(x > 0)})
  ddsnormal <- ddsnormal[filter,]

  
  if (file.exists("sizefactors.1.txt")) {
    sf <- read.table("sizefactors.1.txt", row.names=1, sep="\t", header=TRUE) 
    sizeFactors(ddsnormal) <- c(t(sf))
  } else {
    ddsnormal <- estimateSizeFactors(ddsnormal)
    sizeFactors(ddsnormal) <- sizeFactors(ddsnormal) / mean(sizeFactors(ddsnormal))
    sf <- sizeFactors(ddsnormal) 
    write.table(sf, file="sizefactors.1.txt", quote=FALSE, row.names=TRUE, sep="\t")
  }
 
  # step 2
  if (file.exists("ddsGenes.LRT.RData")) {
    load("ddsGenes.LRT.RData")
  } else {
    ddsnormal <- estimateDispersions(ddsnormal)
    ddsnormal <- nbinomLRT(ddsnormal, reduced = ~ nested.batch, maxit=500)
    if (any(mcols(ddsnormal)$betaConv)) {
      ddsnormal <- ddsnormal[which(mcols(ddsnormal)$betaConv),]
    }
    save(ddsnormal, file="ddsGenes.LRT.RData")
  }
 
  if (file.exists("res.LRT.RData")) {
    load("res.LRT.RData")
  } else {
    res <- results(ddsnormal)
    save(res, file="res.LRT.RData")
  }
 
  FDR = 1.0e-50
  floorPDEG = 1.0e-50
  pval <- res$pvalue
  pval[is.na(pval)] <- 1
  qval <- p.adjust(pval, method = "BH")
  if (sum(qval < FDR) > (floorPDEG * nrow(ddsnormal))) {
      is.DEG <- as.logical(qval < FDR)
  } else {
    is.DEG <- as.logical(rank(pval, ties.method = "min") <= nrow(ddsnormal) * floorPDEG)
  }
 
  # step 3
  nonDEGgenes <- rownames(counts(ddsnormal, normalized=FALSE)[!is.DEG,])
  write.table(nonDEGgenes, file="nonDEGgenes.txt", quote=FALSE, row.names=TRUE, sep="\t")
  ctrlgenesidx <-  match(nonDEGgenes, rownames(cntmatrix))
  ddsGenes <- DESeqDataSetFromMatrix(
    countData = cntmatrix,
    colData = coldata,
    design = ~ nested.batch + type )
  
  ddsGenes <- estimateSizeFactors(ddsGenes, controlGenes = ctrlgenesidx)
 
  #norm.factors <- sizeFactors(ddsGenes)/colSums(cntGenes)
  #norm.factors <- norm.factors/mean(norm.factors)
  #sizeFactors(ddsGenes) <- norm.factors
  write.table(sizeFactors(ddsGenes), file="sizefactors.2.txt", quote=FALSE, row.names=TRUE, sep="\t")
  save(ddsGenes, file="ddsGenes.DEGES.RData")
 
  normalcnts <- counts(ddsGenes, normalized=TRUE)
  idx.nz <- apply(normalcnts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  pdf("multidensity.pdf")
  multidensity( normalcnts[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  pdf("multiecdf.pdf")
  multiecdf( normalcnts[idx.nz ,],
             xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  
 }
 
 
 
 if (file.exists("ddsNew.DEGES.RData")) {
  load(file = "ddsNew.DEGES.RData")
 } else {
 
  cntmatrix = c() 
  if (file.exists("cntmatrix.gene.txt")) {
    cntmatrix = read.table(file="cntmatrix.gene.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
    cntGenes <- cntmatrix
  } else { 
    cntmatrix = matrix(rep(0, nrow(patient_types)*length(genenames)),  nrow=length(genenames))
    for (i in 1:nrow(patient_types)) {
      cnts <- read.table(paste("RawCnts",as.character(patient_types[i,2]),sep="/"), header=TRUE, row.names=1)
      #cntmatrix[,(i*2-1):(i*2)]=as.matrix(cnts)
      cntmatrix[,i]=as.matrix(cnts)
    }
    dim(cntmatrix)
    rownames(cntmatrix) <- genenames
    cntGenes <- cntmatrix
    colnames(cntmatrix) <- coldata$patient
    write.table(cntmatrix, file="cntmatrix.gene.txt", quote=FALSE, row.names=TRUE, sep="\t")
  } 
  
  ddsGenes <- DESeqDataSetFromMatrix(
    countData = cntGenes,
    colData = coldata,
    design = ~ nested.batch + type )
  
  ddsGenes$type <- droplevels(ddsGenes$type)
 
  sf <- read.table("sizefactors.2.txt", row.names=1, sep="\t", header=TRUE)
  sf <- c(t(sf))
  sizeFactors(ddsGenes) <- sf
  
  colnormal = colData(ddsGenes)
  
 # construct counts from discounted elements 
  file.names <- dir("RawCnts", pattern =".discount.ele.cntTable")
  patientIDs <- substr(file.names, 9, 12)
  TEnames <- rownames(read.table(paste("RawCnts",as.character(file.names[1]),sep="/"), header=TRUE, row.names=1))
 
  TEcntmatrix = matrix(rep(0, length(file.names)*length(TEnames)),  nrow=length(TEnames))
  for (i in 1:length(file.names)) {
    cnts <- read.table(paste("RawCnts",as.character(file.names[i]),sep="/"), header=TRUE, row.names=1)
    TEcntmatrix[,i]=as.matrix(cnts)
  }
  dim(TEcntmatrix)
  rownames(TEcntmatrix) <- TEnames
  colnames(TEcntmatrix) <- patientIDs
  
 # merge with existing gene cntmatrix
  t_TEcntmatrix = t(TEcntmatrix)
  t_cntmatrix = t(cntmatrix)
  
  m <- merge(t_cntmatrix, t_TEcntmatrix, by.x="row.names", by.y="row.names")
  new_patients = m[,1]
  new_genenames = colnames(m[,-1])
  new_cntmatrix = t(m[,-1])
  colnames(new_cntmatrix) = new_patients
  write.table(new_cntmatrix, file="cntmatrix.geneTE.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
  stopifnot( colnormal$patient == new_patients ) # need patient ids in same order to use previously estimated sizeFactors
  colnames(new_cntmatrix) = NULL
  ddsnew <- DESeqDataSetFromMatrix(
    countData = new_cntmatrix,
    colData = colnormal,
    design = ~ nested.batch + type )
  ddsnew$type <- droplevels(ddsnew$type)
  # set sizefactors 
  sizeFactors(ddsnew) <- sf
  save(ddsnew, file="ddsNew.DEGES.RData")
  ddsnewcnts <- counts (ddsnew, normalized=TRUE)
  write.table(ddsnewcnts, file="normalizedcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
 
  L1HSnormalized <- ddsnewcnts["L1HS:L1:LINE",]
  write.table(L1HSnormalized, file="L1HS.normalized.txt", quote=FALSE, row.names=TRUE, sep="\t")
  log2colsums <- log2(colSums(ddsnewcnts))
  write.table(log2colsums, file="log2colsums.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
 }
 
 
 
 
 
 if (file.exists("VSTcnts.DEGES.RData")) {
  load(file = "VSTcnts.DEGES.RData")
 } else {
  vsd <- vst(ddsnew, blind = FALSE)
  VSTcnts <- assay(vsd)
  write.table(VSTcnts["L1HS:L1:LINE",], "L1HS.VST.cnts.txt", quote=FALSE, row.names = TRUE)
  save(VSTcnts, file = "VSTcnts.DEGES.RData")
  write.table(VSTcnts, file="VSTcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
  
  genenames = rownames(ddsnew)
  genenames <- gsub("ALR/Alpha", "ALR_Alpha", genenames)
  genenames <- gsub("BSR/Beta", "BSR_Beta", genenames)
  rownames(ddsnew) <- genenames
  coldataNew = colData(ddsnew)
  for (i in 1:nrow(VSTcnts)) {
    fname = paste0("VSTcnts.normal/",genenames[i], ".txt")
    cnts <- matrix(VSTcnts[i,], ncol=1)
    rownames(cnts) <- coldataNew$patient
    colnames(cnts) <- c("patient\tgene")
    write.table(cnts, file=fname, quote=FALSE, sep="\t", row.names=TRUE)
  }
  L1HSnormalized = read.table("L1HS.normalized.txt", row.names=1, sep="\t", check.names = FALSE)
  L1HStable <- cbind.data.frame(colData(ddsnew), L1HSnormalized )
  L1HStable <- cbind.data.frame(L1HStable, VSTcnts["L1HS:L1:LINE",] )
  colnames(L1HStable)[7:8] = c("normalizedcnt", "VSTcnt")
  write.table(L1HStable, file="L1HS.VST.txt", quote=FALSE, row.names=TRUE, sep="\t")
 }
 
 #heatmap
 if (! file.exists("heatmapTop1000Var.pdf")) {
  coldataNew = colData(ddsnew)
  varGenes <- rowVars(VSTcnts)
  topVarianceGenes <- head(order(varGenes, decreasing=T),1000)
  matrix <- VSTcnts[ topVarianceGenes, ]
  matrix <- matrix - rowMeans(matrix)
  # select the 'contrast' you want
  annotation_data <- as.data.frame(coldataNew$type)
  rownames(annotation_data) <- colnames(matrix)
  colnames(annotation_data) <- "type"
  #colors
#colors
  col_vector = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")

  mycolors=list(type = col_vector[1:length(levels(annotation_data$type))])
names(mycolors$type) <- levels(annotation_data$type) 
pdf(file="heatmapTop1000Var.pdf", width=7*10, height=7*12)
par(ps=3)
pheatmap(matrix, 
         annotation_col=annotation_data,
         annotation_colors = mycolors,
         #fontsize = 7
)
  dev.off()
 }
 

  if (file.exists("ddsNew.LRT.RData")) {
    load("ddsNew.LRT.RData")
  } else {
    ddsnew <- estimateDispersions(ddsnew)
    ddsnew <- nbinomLRT(ddsnew, reduced = ~ nested.batch, maxit=500)
    if (any(mcols(ddsnew)$betaConv)) {
      ddsnew <- ddsnew[which(mcols(ddsnew)$betaConv),]
    }
    save(ddsnew, file="ddsNew.LRT.RData")
  }

  #contrasts <- resultsNames(ddsnew)
  coldataNew = colData(ddsnew)
  typelist <- levels(coldataNew$type)
  for (i in 1:length(typelist)) {
    for (j in (i+1):length(typelist)) {
      filename = paste0("resNew.LRT.", i, ".", j, ".RData")
      if (file.exists(filename)) {
        load(filename)
      } else {
        res <- results(ddsnew, contrast=c("type", as.character(typelist[i]), as.character(typelist[j])))
        save(res, file=filename)
      }
    }
  }

 
  for (i in 1:length(typelist)) {
    for (j in (i+1):length(typelist)) {
      resfilename = paste0("resNew.LRT.", i, ".", j, ".RData")
      pdffilename = paste0("MAplot.", typelist[i], ".", typelist[j], ".pdf")
      if (file.exists(resfilename)) {
        load(resfilename)
        pdf(pdffilename)
        plotMA(res)
        dev.off()
      }
    }
  }




# violin plot
L1HS_bytype = read.table("L1HS.VST.txt", header=TRUE, row.names=1, sep="\t")
dodge <- position_dodge(width = 0.6)
ggplot(L1HS_bytype, aes(x=type, y=VSTcnt, fill=condition)) + geom_violin(position = dodge) + geom_boxplot(width=0.1, outlier.colour=NA, position = dodge) + theme_bw()
ggsave('violin.bytype.pdf', dpi=600)



