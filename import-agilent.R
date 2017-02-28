
library(limma)


files <- list.files('AG1/')
files <- files[1:10]
files <- paste0('AG1/',files)

agl_microarray <- read.maimages(files, source="agilent", green.only=TRUE)

######## correction du background

agl_microarray_bc <- backgroundCorrect(agl_microarray,method='subtract')

######### normalization between microarray

agl_microarray_bc_no <- normalizeBetweenArrays(agl_microarray_bc, method="quantile")
par(mfrow=c(2,1))
plotDensities(agl_microarray_bc,from=0,to=16)
plotDensities(agl_microarray_bc_no,from=0,to=16)

mycondition <- c(rep('A',5),rep('B',5))
mydesign <- model.matrix(~mycondition )

fit <- lmFit(agl_microarray_bc_no, design=mydesign)
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
head(resultat.agl,n=20)


resultat.agl[resultat.agl$annotation=='RRM2',]


