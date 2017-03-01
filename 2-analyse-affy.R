library(limma)
library(gcrma)
library(hgu133plus2.db)

files <- list.files('AFFX/')
files <- files[1:10]
files <- paste0('AFFX/',files)

rawfiles <- ReadAffy(filenames=files)
expressionset <- gcrma(rawfiles)
expressionmatrix <- exprs(expressionset)
dim(expressionmatrix)
expressionmatrix[1:10,1:10]
probename <- rownames(expressionset)

columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
annotation <- select(hgu133plus2.db, keys=probename, columns = c("SYMBOL"),keytype='PROBEID')
head(annotation,n=10)
dim(annotation)

idx <- match(probename,annotation$PROBEID)
annotation <- annotation[idx,]
dim(annotation)

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
adj.pvalue.affx <- p.adjust(resultat.affx$pvalue.affx,method='BH')
resultat.affx <- data.frame(resultat.affx,adj.pvalue.affx)

write.table(resultat.affx,'resultat-affx.txt',sep='\t',row.names=F)

head(resultat.affx,n=10)










