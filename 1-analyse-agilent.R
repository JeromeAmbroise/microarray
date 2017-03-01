rm(list=ls())
library(limma)


files <- list.files('1-data/AGL/')
files <- paste0('1-data/AGL/',files)

agl <- read.maimages(files, source="agilent", green.only=TRUE)

######## correction du background

agl_BC <- backgroundCorrect(agl,method='subtract')

######### normalization between microarray

agl_BC_NO <- normalizeBetweenArrays(agl_BC, method="quantile")

pdf('2-result/1a-agilent-normalisation.pdf',width=10,height=10)
par(mfrow=c(2,1))
plotDensities(agl_BC,from=0,to=16)
plotDensities(agl_BC_NO,from=0,to=16)
dev.off()

mycondition <- c(rep('A',5),rep('B',5))
mydesign <- model.matrix(~mycondition )

fit <- lmFit(agl_BC_NO, design=mydesign)
efit <- eBayes(fit)

coefficient.agl <- fit$coefficients
head(coefficient.agl)
coefficient.agl <- coefficient.agl[,2]
pvalue.agl <- efit$p.value
head(pvalue.agl)
pvalue.agl <- pvalue.agl[,2]
annotation <- efit$genes
annotation <- annotation[,'GeneName']

resultat.agl <- data.frame(annotation,coefficient.agl,pvalue.agl)
resultat.agl <- resultat.agl[sort.list(resultat.agl$pvalue.agl),]
adj.pvalue.agl <- p.adjust(resultat.agl$pvalue.agl,method='BH')
resultat.agl <- data.frame(resultat.agl,adj.pvalue.agl)
head(resultat.agl,n=10)

write.table(resultat.agl,'2-result/1-fold-change-agl.txt',sep='\t',row.names=F)


####### graphique


## heatmap

pdf('2-result/1b-agilent-heatmap.pdf',width=8,height=8)
expression <- agl_BC_NO$E[1:100,]
heatmap(expression)
dev.off()

## volcano plot

pdf('2-result/1c-agilent-heatmap.pdf',width=8,height=8)
plot(resultat.agl$coefficient.agl,-log10(resultat.agl$pvalue.agl))
dev.off()

