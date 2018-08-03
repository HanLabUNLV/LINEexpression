library("RColorBrewer")
library("DESeq2")
library("geneplotter")
library("EDASeq")
library("genefilter")
library("pheatmap")

load(file = "ddsNew.DEGES.RData")
load(file = "VSTcnts.DEGES.RData")

lineVSTcnts <- VSTcnts[grep("LINE", rownames(VSTcnts)), ]
LINEcnt <- 2^(lineVSTcnts)
LINEtotal <- colSums(LINEcnt)
L1HScnts <- LINEcnt["L1HS:L1:LINE",]
L1PA2cnts <- LINEcnt["L1PA2:L1:LINE",]
L1PA3cnts <- LINEcnt["L1PA3:L1:LINE",]
LINEyoung <- L1HScnts+L1PA2cnts+L1PA3cnts
LINEold <- LINEtotal - LINEyoung
L1HSVSTcnt <- log2(L1HScnts)
LINEVSTold <- log2(LINEold)

#write.table(LINEVSTold, "LINE.old.VSTcnt.txt", quote=FALSE, row.names=TRUE, sep="\t")
#sineVSTcnts <- VSTcnts[grep("SINE", rownames(VSTcnts)), ]
#hervVSTcnts <- VSTcnts[grep("HERV", rownames(VSTcnts)), ]

#heatmap
coldataNew = colData(ddsnew)
varGenes <- rowVars(lineVSTcnts)
#varGenes <- rowVars(sineVSTcnts)
#varGenes <- rowVars(hervVSTcnts)
topVarianceGenes <- head(order(varGenes, decreasing=T),1000)
matrix <- lineVSTcnts[ topVarianceGenes, ]
#matrix <- sineVSTcnts[ topVarianceGenes, ]
#matrix <- hervVSTcnts[ topVarianceGenes, ]
matrix <- matrix - rowMeans(matrix)
# select the 'contrast' you want
annotation_data <- as.data.frame(coldataNew$type)
rownames(annotation_data) <- colnames(matrix)
colnames(annotation_data) <- "type"

#colors
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")

mycolors=list(type = col_vector[1:length(levels(annotation_data$type))])
names(mycolors$type) <- levels(annotation_data$type) 
#png(filename="LINEheatmapTop1000Var.png", width=6285, height=4860)
#png(filename="SINEheatmapTop1000Var.png", width=6285, height=4860)
#png(filename="HERVheatmapTop1000Var.png", width=6285, height=4860)
pdf(file="LINEheatmapTop1000Var.pdf", width=7*10, height=7*6)
par(ps=3)
pheatmap(matrix,
         annotation_col=annotation_data,
         annotation_colors = mycolors,
         #fontsize = 7
)
dev.off()
