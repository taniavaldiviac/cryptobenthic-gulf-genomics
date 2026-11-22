suppressMessages(library(lme4))
args <- commandArgs(trailingOnly=TRUE)
d <- read.csv(if (length(args)) args[1] else "datos_modelo_v2.csv")
# Expected columns: familia, grupo (C/A cryptobenthic, R conspicuous, O other), ref12S, refMiFish, refteleo,
# refRiaz, ref16S, refCOI, mito (0/1), anio (year of description), occ (GBIF records),
# ms (GBIF MATERIAL_SAMPLE records), ps (PRESERVED_SPECIMEN records), sq (isSequenced records)
# Optional column gbif_rank: rows matched by GBIF above species level get NA for year and record counts
if ("gbif_rank" %in% names(d)) { g <- d$gbif_rank != "SPECIES"; for (v in c("anio","occ","ms","ps","sq")) d[[v]][g] <- NA }
d$grupo3 <- factor(ifelse(d$grupo %in% c("C","A"), "cripto", ifelse(d$grupo=="R","conspicuo","otros")),
                   levels=c("conspicuo","cripto","otros"))
z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
d$anio_s <- (d$anio-1900)/50
d$occ_nodna <- pmax(d$occ - ifelse(is.na(d$ms),0,d$ms), 0)
d$l_occ   <- z(log10(d$occ+1))
d$l_nodna <- z(log10(d$occ_nodna+1))
d$l_ps    <- z(log10(d$ps+1))
elasmo <- c("Myliobatidae","Alopiidae","Rajidae","Scyliorhinidae","Arhynchobatidae","Carcharhinidae","Lamnidae",
            "Cetorhinidae","Narcinidae","Echinorhinidae","Triakidae","Ginglymostomatidae","Gymnuridae","Rhinochimaeridae",
            "Heterodontidae","Hexanchidae","Chimaeridae","Dasyatidae","Odontaspididae","Rhinobatidae","Pristidae",
            "Rhincodontidae","Sphyrnidae","Squatinidae","Torpedinidae","Urotrygonidae","Urolophidae")
d$elasmo <- d$familia %in% elasmo
dd <- d[complete.cases(d[,c("anio_s","l_occ","l_nodna","l_ps")]),]
cat("n =", nrow(dd), "; crypto =", sum(dd$grupo3=="cripto"), "; elasmobranchs =", sum(dd$elasmo), "\n")
cat("Cryptobenthic share of GBIF records that are material samples (median %):",
    round(100*median(dd$ms[dd$grupo3=="cripto"]/pmax(dd$occ[dd$grupo3=="cripto"],1)),1), "\n")
fit <- function(y, rhs, data) glmer(as.formula(paste(y,"~",rhs,"+ (1|familia)")), data, family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))
specs <- list(m0 = "grupo3",
              m1 = "grupo3 + anio_s + l_occ",
              m1_sinADN = "grupo3 + anio_s + l_nodna",
              m1_ejemplares = "grupo3 + anio_s + l_ps")
out <- list()
for (y in c("ref12S","refMiFish","refteleo","refRiaz","ref16S","refCOI","mito")) {
  for (sub in c("todas","sin_elasmo")) {
    dat <- if (sub=="todas") dd else dd[!dd$elasmo,]
    for (nm in names(specs)) {
      if (sub=="sin_elasmo" && nm!="m1_sinADN") next
      m <- tryCatch(fit(y, specs[[nm]], dat), error=function(e) NULL); if (is.null(m)) next
      cf <- summary(m)$coefficients; ci <- confint(m, method="Wald")
      vf <- as.numeric(VarCorr(m)$familia)
      for (term in rownames(cf)[-1]) out[[length(out)+1]] <- data.frame(respuesta=y, datos=sub, modelo=nm, termino=term,
        OR=exp(cf[term,1]), IC95_inf=exp(ci[term,1]), IC95_sup=exp(ci[term,2]), p=cf[term,4],
        ICC_familia=vf/(vf+pi^2/3), AIC=AIC(m), n=nrow(dat))
    }
  }
}
res <- do.call(rbind,out); res[,5:10] <- signif(res[,5:10],3)
write.csv(res, "modelo_sesgo_v2_resultados.csv", row.names=FALSE)
key <- res[res$termino=="grupo3cripto",c("respuesta","datos","modelo","OR","IC95_inf","IC95_sup","p")]
print(key, row.names=FALSE)
eff <- res[grepl("^l_",res$termino) | res$termino=="anio_s", c("respuesta","datos","modelo","termino","OR","IC95_inf","IC95_sup")]
print(eff[eff$datos=="todas",], row.names=FALSE)
