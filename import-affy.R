library(limma)
library(gcrma)

files <- list.files('AFFX/')
files <- files[1:10]
files <- paste0('AFFX/',files)

rawfiles <- ReadAffy(filenames=files)
expressionset <- gcrma(rawfiles)
expressionmatrix <- exprs(expressionset)
expressionmatrix[1:10,1:10]

probename <- rownames(expressionset)


library(hgu133plus2.db)
columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
annotation <- select(hgu133plus2.db, keys=probename, columns = c("SYMBOL"),keytype='PROBEID')

head(annotation,n=10)

idx <- match(probename,annotation$PROBEID)
annotation <- annotation[idx,]

all.equal(annotation$PROBEID,probename)

fData(expressionset) 
fData(expressionset) <- annotation

condition <- c(rep('A',5),rep('B',5))
design <- model.matrix(~condition )
fit <- lmFit(expressionset, design)
efit <- eBayes(fit)

coefficient.affx <- fit$coefficients
head(coefficient.affx)
coefficient.affx <- coefficient.affx[,2]
pvalue.affx <- efit$p.value
head(pvalue.affx)
pvalue.affx <- pvalue.affx[,2]
annotation <- efit$genes
annotation <- annotation[,'SYMBOL']

sum(is.na(annotation))

resultat.affx <- data.frame(annotation,coefficient.affx,pvalue.affx)
dim(resultat.affx)
resultat.affx <- data.frame(na.omit(resultat.affx))
dim(resultat.affx)

resultat.affx <- resultat.affx[sort.list(resultat.affx$pvalue.affx),]
head(resultat.affx,n=20)
resultat.affx[resultat.affx$annotation=='RRM2',]

idx <- match(resultat.agl$annotation,resultat.affx$annotation)

resultat.affx.sort <- resultat.affx[idx,]
compa <- data.frame(resultat.agl,resultat.affx.sort)

plot(compa$coefficient.agl,compa$coefficient.affx)
cor(compa$coefficient.agl,compa$coefficient.affx,use='complete.obs')

cor.test(compa$coefficient.agl,compa$coefficient.affx,use='complete.obs')


