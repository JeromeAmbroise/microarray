idx <- match(resultat.agl$annotation,resultat.affx$annotation)
resultat.affx.sort <- resultat.affx[idx,]
compa <- data.frame(resultat.agl,resultat.affx.sort)

plot(compa$coefficient.agl,compa$coefficient.affx)
cor(compa$coefficient.agl,compa$coefficient.affx,use='complete.obs')