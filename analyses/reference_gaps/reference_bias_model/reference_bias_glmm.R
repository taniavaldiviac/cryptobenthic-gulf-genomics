suppressMessages({library(lme4)})
d <- read.csv("datos_modelo.csv")
# Rows whose GBIF match is genus-level only (column gbif_rank != "SPECIES") get NA for anio and occ
if ("gbif_rank" %in% names(d)) { g <- d$gbif_rank != "SPECIES"; d$anio[g] <- NA; d$occ[g] <- NA }
d$grupo3 <- factor(ifelse(d$grupo %in% c("C","A"), "cripto", ifelse(d$grupo=="R","conspicuo","otros")), levels=c("conspicuo","cripto","otros"))
d$anio_s <- (d$anio-1900)/50          # per 50 years since 1900
d$logocc <- log10(d$occ+1)
d$logocc_s <- (d$logocc-mean(d$logocc,na.rm=TRUE))/sd(d$logocc,na.rm=TRUE)
cat("n =",nrow(d)," complete =",sum(complete.cases(d[,c("anio","occ")]))," families =",length(unique(d$familia)),"\n")
print(aggregate(cbind(ref12S,refCOI,mito,anio,logocc)~grupo3,d,mean))
dd <- d[complete.cases(d[,c("anio_s","logocc_s")]),]
out <- list()
for (y in c("ref12S","refteleo","refMiFish","ref16S","refCOI","mito")) {
  f0 <- as.formula(paste(y,"~ grupo3 + (1|familia)"))
  f1 <- as.formula(paste(y,"~ grupo3 + anio_s + logocc_s + (1|familia)"))
  m0 <- glmer(f0, dd, family=binomial, control=glmerControl(optimizer="bobyqa"))
  m1 <- glmer(f1, dd, family=binomial, control=glmerControl(optimizer="bobyqa"))
  for (nm in c("m0","m1")) { m <- get(nm); cf <- summary(m)$coefficients; ci <- confint(m, method="Wald")
    vf <- as.numeric(VarCorr(m)$familia); icc <- vf/(vf+pi^2/3)
    for (term in rownames(cf)[-1]) out[[length(out)+1]] <- data.frame(respuesta=y, modelo=nm, termino=term,
      OR=exp(cf[term,1]), IC95_inf=exp(ci[term,1]), IC95_sup=exp(ci[term,2]), p=cf[term,4],
      var_familia=vf, ICC_familia=icc, AIC=AIC(m)) }
}
res <- do.call(rbind,out); res[,4:10] <- signif(res[,4:10],3)
print(res, row.names=FALSE)
write.csv(res,"modelo_sesgo_referencia_resultados.csv",row.names=FALSE)
# within-cryptobenthic: which traits predict reference?
dc <- dd[dd$grupo3=="cripto",]
mc <- glmer(ref12S ~ anio_s + logocc_s + (1|familia), dc, family=binomial, control=glmerControl(optimizer="bobyqa"))
cat("\nSolo criptobentonicos (n=",nrow(dc),")\n"); print(round(cbind(OR=exp(fixef(mc)),exp(confint(mc,method="Wald"))[-1,][match(names(fixef(mc)),rownames(confint(mc,method="Wald"))[-1]),],p=summary(mc)$coefficients[,4]),3))
print(VarCorr(mc))
fam <- aggregate(cbind(ref12S,n=1)~familia, dc, sum); fam$prop <- round(fam$ref12S/fam$n,2); print(fam[order(-fam$n),])
