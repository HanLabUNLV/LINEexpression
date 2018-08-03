#library(nlme)
library(MuMIn)
library(lmSupport)
library(qvalue)
library(matrixStats)


generatio <- function(genename)
{
  normalfname <- paste("/home/mirahan/Work/TE-pre/L1HSstats.DEGES/GDC/VSTcnts.normal/", genename, ".txt", sep="")
  
  if (file.info(normalfname)$size == 0) {
    print ("file size zero")
    return (0)
  }
  genen <- read.table(normalfname, header=TRUE)
  #genen <- genen[genen$condition == "normal",]
  if (ncol(genen)==3) genen <- genen[,-2]
  patients <- substr(genen[,1], 1, 4)
  gene <- cbind.data.frame(patients, genen[,2])
  colnames(gene) <- c("patient", "gene_n")
  return (gene)
}


setwd("/home/mirahan/Work/TE-pre/L1HSstats.DEGES/")

L1HS <- read.table("/home/mirahan/Work/TE-pre/L1HSstats.DEGES/L1HS.VST.txt", sep="\t", header=TRUE, row.names=1)
L1HS <- L1HS[L1HS$condition=="normal",]
colnames(L1HS)[7:8] <-c("L1HSnRPM", "L1HS_n")

log2TPMsum <- read.table(file="/home/mirahan/Work/TE-pre/L1HSstats.DEGES/log2colsums.txt", sep="\t", header = TRUE, row.names=1)
log2TPMsum <- cbind(substr(rownames(log2TPMsum), 1, 4), log2TPMsum)
colnames(log2TPMsum) <- c("patient", "log2TPMsum")

oldLINEVST <- read.table(file="/home/mirahan/Work/TE-pre/L1HSstats.DEGES/LINE.old.VSTcnt.txt", sep="\t", header = TRUE, row.names=1)
oldLINEVST <- cbind.data.frame(rownames(oldLINEVST), oldLINEVST)
colnames(oldLINEVST) <- c("patient", "oldLINEVST")

file.names <- dir("/home/mirahan/Work/TE-pre/L1HSstats/GDC/VSTcnts.normal/", pattern =".txt")
file.names.nchar <- sapply(file.names, nchar)
genenames <- substr(file.names, 1, file.names.nchar-4)
genenames <- unique(genenames)

cancertypes = list("BLCA", "BRCA", "CHOL", c("COAD", "READ"), c("ESCA", "STAD"), "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "THCA", "UCEC")
cancerdirs = unlist(lapply(cancertypes, paste, collapse=""))
newdirs = paste("results.oldLINEcov", cancerdirs, sep="/")
sapply(newdirs, dir.create)

THCAradiation = read.table(file="/home/mirahan/Work/TE-pre/L1HSstats.DEGES/THCA_expressionandradiation.txt", sep="\t", header=TRUE)


for (j in 1:length(cancertypes)) {
  cancertype <- cancertypes[[j]]
  cancerdir <- cancerdirs[j]
  pearson_r <- rep(0, length(genenames))
  pearson_pval <- rep(1, length(genenames))
  spearman_rho <- rep(0, length(genenames))
  spearman_pval <- rep(1, length(genenames))
  bestmodelnames <- rep("mod", length(genenames))
  fixed_generatio <- rep(0, length(genenames))
  generatio_coef <- rep(0,length(genenames))
  generatio_pval <- rep(1,length(genenames))
  #R2marginal <- rep(0, length(genenames))
  #R2conditional <- rep(0, length(genenames))
  partialeta2 <- rep(0, length(genenames))
  tested <- rep(0, length(genenames))
  
  for(i in 1:length(genenames)){
    print (i)
    print (genenames[i])
    gene<- generatio(genenames[i])
    if (is.null(dim(gene))) next;

    data <- merge(L1HS, gene, by = "patient")
    data <- merge(data, log2TPMsum, by = "patient")
    data <- merge(data, oldLINEVST, by = "patient")
    if (cancertype == "THCA") {
      data <- merge(data, THCAradiation, by = "patient")
    }
    
    subdata <- data[data$type %in% cancertype,]
    if (nrow(subdata)<8) next;
    if (all(colQuantiles(as.matrix(subdata$gene_n))["25%"]<log2(2))) next;
    subdata$type <- droplevels(subdata$type)
    subdata$batch <- factor(subdata$batch)
    if (length(cancertype)>1) {
      subdata$batch <- droplevels(subdata$batch)
      subdata$nested.batch=as.numeric(subdata$batch)
    }
    subdata$nested.batch = factor(subdata$nested.batch)
    dataratio <- subdata    
    
    test<-cor.test(dataratio$gene_n,dataratio$L1HS_n, method="pearson")
    pearson_r[i] <- test$estimate
    pearson_pval[i] <- test$p.value
    test<-cor.test(dataratio$gene_n,dataratio$L1HS_n, method="spearman")
    spearman_rho[i] <- test$estimate
    spearman_pval[i] <- test$p.value

    AICc = NULL    
    set.seed(0)
    garbage <- rnorm(length(dataratio$L1HS_n))  
    if (cancertype == "THCA") {
      if (length(summary(dataratio$nested.batch))>1) {
        mod0 <- lm(L1HS_n ~ garbage + oldLINEVST, data=dataratio)
        mod1 <- lm(L1HS_n ~ log2TPMsum + oldLINEVST, data=dataratio)
        mod2 <- lm(L1HS_n ~ gene_n + oldLINEVST, data=dataratio)
        mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
        mod4 <- lm(L1HS_n ~ nested.batch + oldLINEVST, data=dataratio)
        mod5 <- lm(L1HS_n ~ log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
        mod6 <- lm(L1HS_n ~ gene_n + nested.batch + oldLINEVST, data=dataratio)
        mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
        mod8 <- lm(L1HS_n ~ radiation + oldLINEVST, data=dataratio)
        mod9 <- lm(L1HS_n ~ log2TPMsum + radiation + oldLINEVST, data=dataratio)
        mod10 <- lm(L1HS_n ~ gene_n + radiation + oldLINEVST, data=dataratio)
        mod11 <- lm(L1HS_n ~ gene_n + log2TPMsum + radiation + oldLINEVST, data=dataratio)
        mod12 <- lm(L1HS_n ~ nested.batch + radiation + oldLINEVST, data=dataratio)
        mod13 <- lm(L1HS_n ~ log2TPMsum + nested.batch + radiation + oldLINEVST, data=dataratio)
        mod14 <- lm(L1HS_n ~ gene_n + nested.batch + radiation + oldLINEVST, data=dataratio)
        mod15 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch + radiation + oldLINEVST, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15)
      } else {
        mod0 <- lm(L1HS_n ~ garbage + oldLINEVST, data=dataratio)
        mod1 <- lm(L1HS_n ~ log2TPMsum + oldLINEVST, data=dataratio)
        mod2 <- lm(L1HS_n ~ gene_n + oldLINEVST, data=dataratio)
        mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
        mod4 <- lm(L1HS_n ~ radiation + oldLINEVST, data=dataratio)
        mod5 <- lm(L1HS_n ~ log2TPMsum + radiation + oldLINEVST, data=dataratio)
        mod6 <- lm(L1HS_n ~ gene_n + radiation + oldLINEVST, data=dataratio)
        mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + radiation + oldLINEVST, data=dataratio)
        AICc<-model.sel(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      }      
    } else {
    if (length(summary(dataratio$nested.batch))>1) {
      mod0 <- lm(L1HS_n ~ garbage + oldLINEVST, data=dataratio)
      mod1 <- lm(L1HS_n ~ log2TPMsum + oldLINEVST, data=dataratio)
      mod2 <- lm(L1HS_n ~ gene_n + oldLINEVST, data=dataratio)
      mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
      mod4 <- lm(L1HS_n ~ nested.batch + oldLINEVST, data=dataratio)
      mod5 <- lm(L1HS_n ~ log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
      mod6 <- lm(L1HS_n ~ gene_n + nested.batch + oldLINEVST, data=dataratio)
      mod7 <- lm(L1HS_n ~ gene_n + log2TPMsum + nested.batch + oldLINEVST, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3,mod4, mod5, mod6, mod7)
    } else {
      mod0 <- lm(L1HS_n ~ garbage + oldLINEVST, data=dataratio)
      mod1 <- lm(L1HS_n ~ log2TPMsum + oldLINEVST, data=dataratio)
      mod2 <- lm(L1HS_n ~ gene_n + oldLINEVST, data=dataratio)
      mod3 <- lm(L1HS_n ~ gene_n + log2TPMsum + oldLINEVST, data=dataratio)
      AICc<-model.sel(mod0,mod1,mod2,mod3)
    }
    }
    
    bestmodel <- eval(getCall(AICc, 1))
    bestmodelnames[i] <- rownames(AICc)[1]
    sm <- summary(bestmodel)
    fixed_generatio[i] <- any(rownames(sm$coefficients)=="gene_n")
    if (fixed_generatio[i]) {
      plot(bestmodel$model[,"gene_n"], resid(bestmodel), ylab="residuals", xlab="gene_n")
      generatio_coef[i] <- sm$coefficients["gene_n","Estimate"]
      generatio_pval[i] <- sm$coefficients["gene_n","Pr(>|t|)"]
      if (genenames[i] != "L1HS:L1:LINE" && !any(sm$aliased) && deviance(bestmodel) > 1.0e-08) {
        effectsize <- modelEffectSizes(bestmodel)
        partialeta2[i] <- effectsize$Effects["gene_n","pEta-sqr"]
      } else {
        partialeta2[i] <- NA
      }
    }
    tested[i] <- 1
    
  }

  testedidx = which(tested==1)
  if (length(testedidx) > 0) {
  pearson_qval = qvalue(pearson_pval[testedidx])$qvalues
  spearman_qval = qvalue(spearman_pval[testedidx])$qvalues
  filename = paste("results.oldLINEcov", cancerdir, "results_cor.txt", sep="/")
  results <- cbind.data.frame(genenames[testedidx], pearson_r[testedidx], pearson_pval[testedidx], pearson_qval, spearman_rho[testedidx], spearman_pval[testedidx], spearman_qval)
  colnames(results) <- c("genenames", "pearson_r" , "pearson_pval", "pearson_qval" ,"spearman_rho", "spearman_pval", "spearman_qval" )
  results <- results[order(-results[,"spearman_rho"]),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  generatio_qval = qvalue(generatio_pval[testedidx])$qvalues 
  filename = paste("results.oldLINEcov", cancerdir, "results_lm.txt", sep="/")
  results <- cbind.data.frame(genenames[testedidx], fixed_generatio[testedidx], bestmodelnames[testedidx], generatio_coef[testedidx], generatio_pval[testedidx], generatio_qval, partialeta2[testedidx])
  colnames(results) <- c("genenames", "fixed_generatio", "bestmodelnames", "generatio_coef", "generatio_pval", "generatio_qval", "partialeta2")
  results <- results[order(-results$partialeta2),]
  write.table(results, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  
  
  filename = paste("results.oldLINEcov", cancerdir, "results_positive_lm.txt", sep="/")
  positive<-as.data.frame(results[results$fixed_generatio==1,])
  positive[,4] <- as.numeric(as.character(positive[,4]))
  positive[,5] <- as.numeric(as.character(positive[,5]))
  positive[,6] <- as.numeric(as.character(positive[,6]))
  positive <- positive[order(-positive$partialeta2),]
  write.table(positive, file=filename, sep="\t", row.names = FALSE, quote=FALSE)
  }
}




coexpressed_plus = c()
coexpressed_minus = c()
for (j in 1:length(cancertypes)) {
  type <- cancertypes[j]
  cancerdir <- cancerdirs[j]
  fname = paste("results.oldLINEcov",cancerdir, "results_positive_lm.txt", sep="/")
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
}
write.table(coexpressed_plus, file="results.oldLINEcov/coexpressed_plus.txt", sep="\t", quote = FALSE, row.names = FALSE )
write.table(coexpressed_minus, file="results.oldLINEcov/coexpressed_minus.txt", sep="\t", quote = FALSE, row.names = FALSE )



