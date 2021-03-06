
library("RDAVIDWebService")
library(org.Hs.eg.db)
library(annotate)


getClusterSymbols <- function (termCluster) {
  N= nrow(summary(termCluster))
  clustersymbols <- vector("list", N)
  for (i in 1:N ) {
    clusterTerms<-members(termCluster)[[i]]
    genes <- unique(unlist(strsplit(clusterTerms$Genes, ", ")))
    symbols <- getSYMBOL(genes, data='org.Hs.eg')
    clustersymbols[[i]] <- symbols
  }
  return (clustersymbols)
}


runDavid <- function(resultdir) {
  david<-DAVIDWebService(email="mira.han@unlv.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  setAnnotationCategories(david, getDefaultCategoryNames(david))
  minuslist <- read.table(paste(resultdir,"duplicate.minus.rsquared",sep="/"), header=FALSE)
  result<-addList(david, minuslist[,2], idType="ENTREZ_GENE_ID", listName="minusduplicates", listType="Gene")
  termClusterMinus<-getClusterReport(david, type="Term")
  symbolfile = "termClusterSymbols.minus.tab"
  if (file.exists(symbolfile)) file.remove(symbolfile)
  if (length(members(termClusterMinus)) > 0) {
    symbolsMinus <- getClusterSymbols(termClusterMinus)
    symbolsMinus <- lapply(symbolsMinus, sort)
    lapply(symbolsMinus, write, paste(resultdir,symbolfile, sep="/"), append=TRUE, ncolumns=1000)
    getClusterReportFile(david, type="Term", fileName=paste(resultdir,"termClusterReport.minus.tab",sep="/"))
    getFunctionalAnnotationChartFile(david, fileName=paste(resultdir,"FAChart.minus.tab",sep="/"))
  }

  #davidGODag<-DAVIDGODag(members(termClusterMinus)[[1]],pvalueCutoff=0.1, "CC")
  #plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")

  pluslist <- read.table(paste(resultdir,"duplicate.plus.rsquared",sep="/"), header=FALSE)
  result<-addList(david, pluslist[,2], idType="ENTREZ_GENE_ID", listName="plusduplicates", listType="Gene")
  termClusterPlus<-getClusterReport(david, type="Term")
  symbolfile = "termClusterSymbols.plus.tab"
  if (file.exists(symbolfile)) file.remove(symbolfile)
  if (length(members(termClusterPlus)) > 0) { 
    symbolsPlus <- getClusterSymbols(termClusterPlus)
    symbolsPlus <- lapply(symbolsPlus, sort)
    lapply(symbolsPlus, write, paste(resultdir,symbolfile, sep="/"), append=TRUE, ncolumns=1000)
    getClusterReportFile(david, type="Term", fileName=paste(resultdir,"termClusterReport.plus.tab",sep="/"))
    getFunctionalAnnotationChartFile(david, fileName=paste(resultdir,"FAChart.plus.tab",sep="/"))
  }

}

setwd("/home/mirahan/Work/TE-pre/L1HSstats.fewertissues")

resultdir = "results.oldLINEcov"
runDavid(resultdir)
resultdir = "results.oldLINE"
runDavid(resultdir)
resultdir = "results.VSTnormal"
runDavid(resultdir)

