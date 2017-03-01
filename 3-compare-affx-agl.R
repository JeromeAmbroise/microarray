
resultat.agl <- read.table('2-result/1-fold-change-agl.txt',sep='\t',header=T)
resultat.affx <- read.table('2-result/2-fold-change-affx.txt',sep='\t',header=T)

head(resultat.agl)
head(resultat.affx)

idx <- match(resultat.agl$annotation,resultat.affx$annotation)
resultat.affx <- resultat.affx[idx,]

head(resultat.agl)
head(resultat.affx)

comparison <- data.frame(resultat.agl,resultat.affx)

plot(comparison$coefficient.agl,comparison$coefficient.affx)
cor(comparison$coefficient.agl,comparison$coefficient.affx,use='complete.obs')