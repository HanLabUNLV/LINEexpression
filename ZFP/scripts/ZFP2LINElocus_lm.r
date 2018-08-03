#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)
library(Xmisc)

args = commandArgs(trailingOnly = TRUE)

#Enter the L1family, ZFPgene, and cancertype arguments exactly how they appear on the L1 filename/ZFPgene filename/cancertype column.
#The directories are hard coded, so they may have to be modified if you're running different tests.

#ZNF662|389114_L2a:L2:LINE_KIRP_coexpressed_minus.txt
L1family = "L2a:L2:LINE"
ZFPgene = "ZNF662|389114"
cancertype = "KIRP"
if (length(args) > 0) {
L1family = args[1]
ZFPgene = args[2]
cancertype = args[3]
}

cancerdir = cancertype
cancername = cancertype

if (cancertype == "CHOLLIHC")
{
  cancertype <- c("CHOL", "LIHC")
}
if (cancertype == "COADREAD")
{
  cancertype <- c("COAD", "READ")
}
if (cancertype == "ESCASTAD")
{
  cancertype <- c("ESCA", "STAD")
}
if (cancertype == "KICHKIRCKIRP")
{
  cancertype <- c("KICH", "KIRC", "KIRP")
}
if (cancertype == "LUADLUSC")
{
  cancertype <- c("LUAD", "LUSC")
}

generatio <- function(genename)
{
  normalfname <- paste("/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/ZincFingerGenes/", genename, ".txt", sep="")

  if (file.info(normalfname)$size == 0) {
    print ("file size zero")
    return (0)
  }
  gene <- read.table(normalfname, header=TRUE)
  colnames(gene) <- c("patient", "gene_n")

  gene$patient <- rstrip(gene$patient, char="_normal")
  return (gene)
}

readL1 <- function(L1dup, p2ct2batch)
{
  L1HS <- read.table(paste("/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/LINEInstances/", L1family, "/", L1dup, ".txt", sep=""), header=TRUE)
  colnames(L1HS)[2] <-"L1HS_n"
  L1HS <- merge(L1HS, p2ct2batch, by="patient")
  return (L1HS)
}

setwd("/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/")

#L1familyname <- rstrip(basename(L1file), char = ".txt")


 batchnormal<- read.table("BatchData.normal.tsv", header = TRUE, sep="\t")
 barcode <- strsplit(as.character(batchnormal[,1]), "-")
 barcode <- matrix(unlist(barcode), ncol = 7, byrow = TRUE)
 p2batchnormal <- cbind.data.frame(barcode[,3], batchnormal[,2])
 colnames(p2batchnormal) <- c("patient", "batch")


 patient2cancer <- read.table("patientToCancerMappings.tsv")
 barcode <- strsplit(as.character(patient2cancer[,1]), "-")
 barcode <- matrix(unlist(barcode), ncol = 3, byrow = TRUE)
 p2ct <- cbind.data.frame(barcode[,3], patient2cancer[,2])
 colnames(p2ct) <- c("patient", "type")

 p2ct2batch <- merge(p2ct, p2batchnormal, by="patient")


log2TPMsum <- read.table(file="/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/log2colsums.txt", sep="\t", header = TRUE)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")
log2TPMsum$patient <- rstrip(log2TPMsum$patient, char="_normal")

#instead of using a list of genes, we will use a list of L1 dups and one gene
file.names <- dir(paste("/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/LINEInstances/", L1family, sep=""), pattern =".txt")
#file.names.nchar <- sapply(file.names, nchar)
L1names <- rstrip(file.names, char=".txt")
L1names <- unique(L1names)
#cancertypes = list("BLCA", "BRCA", "CHOL", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "THCA", "UCEC")
cancertypes = list("BLCA", "BRCA", c("CHOL", "LIHC"), c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", c("KICH", "KIRC", "KIRP"), c("LUAD", "LUSC"), "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
newdirs = paste("InstanceResults", cancerdirs, sep="/")
sapply(newdirs, dir.create)

THCAradiation = read.table(file="/DataDrives/dd4/mirahan/TE-pre/L1HSstats.fewertissues/ZFP/THCA_expressionandradiation.txt", sep="\t", header=TRUE)


#for (j in 1:length(cancertypes)) {
#cancertype <- cancertypes[[j]]
pearson_r <- rep(0, length(L1names))
pearson_pval <- rep(1, length(L1names))
#spearman_rho <- rep(0, length(L1names))
#spearman_pval <- rep(1, length(L1names))
bestmodelnames <- rep("mod", length(L1names))
fixed_generatio <- rep(0, length(L1names))
generatio_coef <- rep(0,length(L1names))
generatio_pval <- rep(1,length(L1names))
#R2marginal <- rep(0, length(L1names))
#R2conditional <- rep(0, length(L1names))
partialeta2 <- rep(0, length(L1names))
tested <- rep(0, length(L1names))
for(i in 1:length(L1names)){
  print (i)
  print (L1names[i])
  print (cancertype)
  gene<- generatio(ZFPgene)
  L1HS<-readL1(L1names[i], p2ct2batch)
  if (is.null(dim(gene))) {
    print("null dimensions for gene")
    next;
  }
  data <- merge(L1HS, gene, by = "patient")
  data <- merge(data, log2TPMsum, by = "patient")
  if (cancertype == "THCA") {
    data <- merge(data, THCAradiation, by = "patient")
  }

  subdata <- data[data$type %in% cancertype,]

  if (all(sapply(subdata$L1HS_n, identical, subdata$L1HS_n[1]))) {
    print("identical L1HS count")
    next;
  }

  if (nrow(subdata)<8) {
    print("less than 8 rows")
    next;
  }
  if (all(colQuantiles(as.matrix(subdata$gene_n))["25%"]<log2(2))) {
    print("col quantiles")
    next;
  }


  subdata$type <- droplevels(subdata$type)
  subdata$batch <- factor(subdata$batch)
  subdata$batch <- droplevels(subdata$batch)
  subdata <- cbind.data.frame(subdata, nested.batch=as.numeric(subdata$batch))
  subdata$nested.batch = factor(subdata$nested.batch)
  dataratio <- subdata

  test<-cor.test(dataratio$gene_n,dataratio$L1HS_n, method="pearson")
  pearson_r[i] <- test$estimate
  pearson_pval[i] <- test$p.value
  #test<-cor.test(dataratio$gene_n,dataratio$L1HS_n, method="spearman")
  #spearman_rho[i] <- test$estimate
  #spearman_pval[i] <- test$p.value

  AICc = NULL
  set.seed(0)
  garbage <- rnorm(length(dataratio$L1HS_n))
  if (cancertype == "THCA") {
    if (length(summary(dataratio$nested.batch))>1) {
      mod0 <- lm(L1HS_n ~ garbage, data=dataratio)
      mod1 <- lm(L1HS_n ~ log2TPMsum, data=dataratio)
      mod2 <- lm(L1HS_n ~ gene_n, data=dataratio)
      mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum, data=dataratio)
      mod4 <- lm(L1HS_n ~ nested.batch, data=dataratio)
      mod5 <- lm(L1HS_n ~ log2TPMsum + nested.batch, data=dataratio)
      mod6 <- lm(L1HS_n ~ gene_n + nested.batch, data=dataratio)
      mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
      mod8 <- lm(L1HS_n ~ radiation, data=dataratio)
      mod9 <- lm(L1HS_n ~ log2TPMsum + radiation, data=dataratio)
      mod10 <- lm(L1HS_n ~ gene_n + radiation, data=dataratio)
      mod11 <- lm(L1HS_n ~ gene_n + log2TPMsum + radiation, data=dataratio)
      mod12 <- lm(L1HS_n ~ nested.batch + radiation, data=dataratio)
      mod13 <- lm(L1HS_n ~ log2TPMsum + nested.batch + radiation, data=dataratio)
      mod14 <- lm(L1HS_n ~ gene_n + nested.batch + radiation, data=dataratio)
      mod15 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch + radiation, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
    } else {
      mod0 <- lm(L1HS_n ~ garbage, data=dataratio)
      mod1 <- lm(L1HS_n ~ log2TPMsum, data=dataratio)
      mod2 <- lm(L1HS_n ~ gene_n, data=dataratio)
      mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum, data=dataratio)
      mod4 <- lm(L1HS_n ~ radiation, data=dataratio)
      mod5 <- lm(L1HS_n ~ log2TPMsum + radiation, data=dataratio)
      mod6 <- lm(L1HS_n ~ gene_n + radiation, data=dataratio)
      mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + radiation, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
    }
  } else {
  if (length(summary(dataratio$nested.batch))>1) {
    mod0 <- lm(L1HS_n ~ garbage, data=dataratio)
    mod1 <- lm(L1HS_n ~ log2TPMsum, data=dataratio)
    mod2 <- lm(L1HS_n ~ gene_n, data=dataratio)
    mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum, data=dataratio)
    mod4 <- lm(L1HS_n ~ nested.batch, data=dataratio)
    mod5 <- lm(L1HS_n ~ log2TPMsum + nested.batch, data=dataratio)
    mod6 <- lm(L1HS_n ~ gene_n + nested.batch, data=dataratio)
    mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch, data=dataratio)
    AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
  } else {
    mod0 <- lm(L1HS_n ~ garbage, data=dataratio)
    mod1 <- lm(L1HS_n ~ log2TPMsum, data=dataratio)
    mod2 <- lm(L1HS_n ~ gene_n, data=dataratio)
    mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum, data=dataratio)
    AICc<-model.sel(mod0,mod1,mod2,mod3)
  }
  }

  bestmodel <- eval(getCall(AICc, 1))
  bestmodelnames[i] <- rownames(AICc)[1]
  sm <- summary(bestmodel)
  fixed_generatio[i] <- any(rownames(sm$coefficients)=="gene_n")
  if (fixed_generatio[i]) {
    #plotfile = paste("/DataDrives/shared_107_dd4/ZincFingers/InstanceResults/Plots/", ZFPgene, "_", L1family, "_", cancername, "_", "plots.pdf", sep="")
    #pdf(file=plotfile)
    #plot(bestmodel$model[,"gene_n"], resid(bestmodel), ylab="residuals", xlab="gene_n")
    #dev.off()
    generatio_coef[i] <- sm$coefficients["gene_n","Estimate"]
    generatio_pval[i] <- sm$coefficients["gene_n","Pr(>|t|)"]
    if (!any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
      effectsize <- modelEffectSizes(bestmodel)
      partialeta2[i] <- effectsize$Effects["gene_n","pEta-sqr"]
    } else {
      generatio_coef[i] <- 0
      generatio_pval[i] <- 1
      partialeta2[i] <- NA
      fixed_generatio[i] <- 0
    }
  }
  tested[i] <- 1

}

testedidx = which(tested==1)
if (length(testedidx) > 0) {
pearson_qval = p.adjust(pearson_pval[testedidx])
#spearman_qval = p.adjust(spearman_pval[testedidx])
filename = paste("InstanceResults/", cancerdir, "/", ZFPgene, "_", L1family, "_results_cor.txt", sep="")
results <- cbind.data.frame(L1names[testedidx], pearson_r[testedidx], pearson_pval[testedidx], pearson_qval)
colnames(results) <- c("L1names", "pearson_r" , "pearson_pval", "pearson_qval" )
#results <- results[order(-results[,"spearman_rho"]),]
write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)

generatio_qval = p.adjust(generatio_pval[testedidx])
filename = paste("InstanceResults/", cancerdir, "/", ZFPgene, "_", L1family, "_results_lm.txt", sep="")
results <- cbind.data.frame(L1names[testedidx], fixed_generatio[testedidx], bestmodelnames[testedidx], generatio_coef[testedidx], generatio_pval[testedidx], generatio_qval, partialeta2[testedidx])
colnames(results) <- c("L1names", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "generatio_qval", "partialeta2")
results <- results[order(-results$partialeta2),]
write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)


filename = paste("InstanceResults/", cancerdir, "/", ZFPgene, "_", L1family, "_results_positive_lm.txt", sep="")
positive<-as.data.frame(results[results$fixed_generatio==1,])
positive[,4] <- as.numeric(as.character(positive[,4]))
positive[,5] <- as.numeric(as.character(positive[,5]))
positive[,6] <- as.numeric(as.character(positive[,6]))
positive <- positive[order(-positive$partialeta2),]
write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
}
#}




coexpressed_plus = c()
coexpressed_minus = c()
#for (j in 1:length(cancertypes)) {
fname = paste("InstanceResults/", cancerdir, "/", ZFPgene, "_", L1family, "_results_positive_lm.txt", sep="")
if (file.exists(fname)) {
  eval <- read.table(fname, sep="\t", header=TRUE)
  eval <- eval[order(-eval$partialeta2),]
  minuscoef <- eval[eval$generatio_coef<0,]
  pluscoef <- eval[eval$generatio_coef>0,]
  dim(minuscoef)
  dim(pluscoef)
  sigplus <- pluscoef[pluscoef$generatio_qval<0.0001,]
  sigminus <- minuscoef[minuscoef$generatio_qval<0.0001,]
  coexpressed_plus = rbind.data.frame(coexpressed_plus, cbind.data.frame(rep(cancerdir, nrow(sigplus)), sigplus))
  coexpressed_minus = rbind.data.frame(coexpressed_minus, cbind.data.frame(rep(cancerdir, nrow(sigminus)), sigminus))
}
#}
filename = paste("InstanceResults/Coexpressed/", ZFPgene, "_", L1family, "_", cancername, "_coexpressed_plus.txt", sep="")
write.table(coexpressed_plus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
filename = paste("InstanceResults/Coexpressed/", ZFPgene, "_", L1family, "_", cancername, "_coexpressed_minus.txt", sep="")
write.table(coexpressed_minus, file=filename, sep="\t", quote = FALSE, row.names = FALSE )
