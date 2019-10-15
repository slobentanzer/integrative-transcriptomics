rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#GSE5392

library(BiocParallel)
registered()
library(lme4)
library(GEOquery)
library(Biobase)
library(BiocParallel)
library(affy)
library(hugene10stv1cdf)
library(limma)


if(!file.exists("./working_data/subset/bahn_raw_eset.rds")) {
  
  #META####
  if (!file.exists("05_subset_analysis/06_bahn/gpl96.rds"))
  {
    gpl96 <- getGEO("GSE5392")
    gpl96 <- gpl96[[1]]
    saveRDS(gpl96, file = "05_subset_analysis/06_bahn/gpl96.rds")
  } else
    gpl96 <- readRDS(file = "05_subset_analysis/06_bahn/gpl96.rds")
  
  #format gpl96####
  pheno <- pData(gpl96)
  head(pheno)
  disease.state <- pheno$characteristics_ch1
  disease.state <- strsplit(as.character(disease.state), ":", fixed = T)
  disease.state <- as.factor(unlist(lapply(disease.state, "[", 2)))
  levels(disease.state) <- c("BD", "CTL")
  tissue <- pheno$title
  tissue <- strsplit(as.character(tissue), "_", fixed = T)
  tissue <- as.factor(unlist(lapply(tissue, "[", 2)))
  hemi <- pheno$characteristics_ch1.7
  hemi <- strsplit(as.character(hemi), ":", fixed = T)
  hemi <- as.factor(unlist(lapply(hemi, "[", 2)))
  sex <- pheno$characteristics_ch1.2
  sex <- strsplit(as.character(sex), ":", fixed = T)
  sex <- as.factor(unlist(lapply(sex, "[", 2)))
  age <- pheno$characteristics_ch1.1
  age <- strsplit(as.character(age), ":", fixed = T)
  age <- as.numeric(unlist(lapply(age, "[", 2)))
  onset <- pheno$characteristics_ch1.3
  onset <- strsplit(as.character(onset), ":", fixed = T)
  onset <- as.numeric(unlist(lapply(onset, "[", 2)))
  onset[disease.state == "CTL"] <- 0
  duration <- pheno$characteristics_ch1.4
  duration <- strsplit(as.character(duration), ":", fixed = T)
  duration <- as.numeric(unlist(lapply(duration, "[", 2)))
  duration[disease.state == "CTL"] <- 0
  ph <- pheno$characteristics_ch1.5
  ph <- strsplit(as.character(ph), ":", fixed = T)
  ph <- as.numeric(unlist(lapply(ph, "[", 2)))
  pmi <- pheno$characteristics_ch1.6
  pmi <- strsplit(as.character(pmi), ":", fixed = T)
  pmi <- as.numeric(unlist(lapply(pmi, "[", 2)))
  fluphenazine.mg.equiv <- pheno$characteristics_ch1.8
  fluphenazine.mg.equiv <- strsplit(as.character(fluphenazine.mg.equiv), ":", fixed = T)
  fluphenazine.mg.equiv <- as.numeric(unlist(lapply(fluphenazine.mg.equiv, "[", 2)))
  fluphenazine.mg.equiv[disease.state == "CTL"] <- 0
  fluphenazine.mg.equiv[is.na(fluphenazine.mg.equiv)] <- 0
  lithium <- pheno$characteristics_ch1.9
  lithium <- strsplit(as.character(lithium), ":", fixed = T)
  lithium <- as.factor(unlist(lapply(lithium, "[", 2)))
  valproate <- pheno$characteristics_ch1.10
  valproate <- strsplit(as.character(valproate), ":", fixed = T)
  valproate <- as.factor(unlist(lapply(valproate, "[", 2)))
  electrotherapy <- pheno$characteristics_ch1.11
  electrotherapy <- strsplit(as.character(electrotherapy), ":", fixed = T)
  electrotherapy <- as.factor(unlist(lapply(electrotherapy, "[", 2)))
  suicide <- pheno$characteristics_ch1.12
  suicide <- strsplit(as.character(suicide), ":", fixed = T)
  suicide <- as.factor(unlist(lapply(suicide, "[", 2)))
  drug.abuse.0.to.6 <- pheno$characteristics_ch1.13
  drug.abuse.0.to.6 <- strsplit(as.character(drug.abuse.0.to.6), ":", fixed = T)
  drug.abuse.0.to.6 <- as.numeric(unlist(lapply(drug.abuse.0.to.6, "[", 3)))
  drug.abuse.0.to.6[is.na(drug.abuse.0.to.6)] <- 0
  alcohol.abuse.0.to.6 <- pheno$characteristics_ch1.14
  alcohol.abuse.0.to.6 <- strsplit(as.character(alcohol.abuse.0.to.6), ":", fixed = T)
  alcohol.abuse.0.to.6 <- as.numeric(unlist(lapply(alcohol.abuse.0.to.6, "[", 3)))
  alcohol.abuse.0.to.6[is.na(alcohol.abuse.0.to.6)] <- 0
  
  meta <- data.frame(disease.state, tissue, age, sex, hemi, onset, duration, fluphenazine.mg.equiv, lithium, valproate, electrotherapy, suicide,
                     drug.abuse.0.to.6, alcohol.abuse.0.to.6, ph, pmi)
  rownames(meta) <- rownames(pheno)
  
  #READ AFFY####
  library(affy)
  affy.folder <- "./raw_data/datasets/Bahn_GSE5392/"
  affy.files <- paste0(affy.folder, list.files(affy.folder))
  affy.files <- affy.files[grep("cel", affy.files)]
  raw_read <- ReadAffy(filenames = affy.files)
  eset <- rma(raw_read)
  
  eset_exprs <- exprs(eset)
  names <- colnames(eset_exprs)
  names <- unlist(lapply(strsplit(names, ".", fixed = T), "[", 1))
  colnames(eset_exprs) <- names
  plotMA(eset_exprs)
  
  idx <- match(colnames(eset_exprs), rownames(meta))
  meta = meta[idx,]
  
  RNAdeg = AffyRNAdeg(raw_read)
  plotAffyRNAdeg(RNAdeg)
  meta$rna.deg = RNAdeg$slope
  
  batch <- factor(substring(protocolData(raw_read)$ScanDate,1,8))
  batch <- data.frame(names, batch)
  
  idx <- match(rownames(meta), batch$names)
  meta$batch <- batch$batch[idx]
  meta$batch <- factor(meta$batch)
  
  mds <- cmdscale(dist(t(eset_exprs)))
  plot(mds, col=as.numeric(as.factor(meta$disease.state)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$disease.state), col=c(1:2), pch=16, cex=0.8)
  
  plot(mds, col=as.numeric(as.factor(meta$tissue)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$tissue), col=c(1:2), pch=16, cex=0.8)
  
  plot(mds, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$batch), col=seq(length(levels(meta$batch))), pch=16, cex=0.8)
  
  plot(mds, col=as.numeric(as.factor(meta$sex)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$sex), col=seq(length(levels(meta$sex))), pch=16, cex=0.8)
  
  plot(mds, col=as.numeric(as.factor(meta$hemi)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$hemi), col=seq(length(levels(meta$hemi))), pch=16, cex=0.8)
  
  rownames(meta) == colnames(eset_exprs)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  eset
  
  saveRDS(eset, "./working_data/subset/bahn_raw_eset.rds")
}
#QC, correlation####

if(!file.exists("./working_data/subset/bahn_regressed_combat.rds")) {
  eset <- readRDS("./working_data/subset/bahn_raw_eset.rds")
  meta <- pData(eset)
  
  #batches
  table(meta$batch) == 1
  idx <-
    which(meta$batch %in% names(table(meta$batch))[table(meta$batch) == 1])
  meta[idx, ]
  eset <- eset[, -idx]
  meta <- pData(eset)
  meta$batch <- factor(meta$batch)
  
  #groups
  plot(meta$disease.state)
  head(meta)
  
  #colinearities
  for (i in 2:length(meta)) {
    # print(i)
    cor <- anova(lm(as.numeric(meta[, i]) ~ meta$disease.state))
    plot(meta[, i] ~ meta$disease.state,
         main = paste(colnames(meta)[i], "=", cor$`Pr(>F)`[1]))
    if (cor$`Pr(>F)`[1] < .05)
      print(colnames(meta)[i])
  }
  # ph, pmi
  
  #CORRECTION####
  library(sva)
  meta <- pData(eset)
  meta$batch <- factor(meta$batch)
  meta$disease.state <- factor(meta$disease.state)
  plot(meta$disease.state ~ meta$batch)
  model = model.matrix( ~ disease.state, data = meta)
  eset_exprs = ComBat(exprs(eset), batch = meta$batch, mod = model)
  
  par(mfrow = c(1, 2))
  mds <- cmdscale(dist(t(exprs(eset))))
  plot(mds,
       col = as.numeric(as.factor(meta$batch)),
       pch = 16,
       main = "MDS Plot")
  legend(
    "topleft",
    levels(meta$batch),
    col = c(1:3),
    pch = 16,
    cex = 0.8
  )
  
  mds.combat <- cmdscale(dist(t(eset_exprs)))
  plot(
    mds.combat,
    col = as.numeric(as.factor(meta$batch)),
    pch = 16,
    main = "MDS Plot ComBat"
  )
  legend(
    "topleft",
    levels(meta$batch),
    col = c(1:3),
    pch = 16,
    cex = 0.8
  )
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  #regression, ph, pmi
  
  eset_exprs <- exprs(eset)
  X <-
    model.matrix( ~ disease.state + tissue + hemi + sex + age + pmi + ph +
                    rna.deg, data = meta)
  Y <- eset_exprs
  beta = (solve(t(X) %*% X) %*% t(X)) %*% t(Y)
  
  reg_corr <- (as.matrix(X[, 6:9]) %*% (as.matrix(beta[6:9, ])))
  eset_exprs <- eset_exprs - t(reg_corr)
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  saveRDS(eset, "./working_data/subset/bahn_regressed_combat.rds")
} 


1
if(F){
eset <- readRDS("./working_data/subset/bahn_regressed_combat.rds")
meta <- pData(eset)

#MDS
par(mfrow = c(1,1))
mds <- cmdscale(dist(t(exprs(eset))), eig = T)
plot(mds$points)
pc1 <- mds$eig[1]^2 / sum(mds$eig^2)
pc2 <- mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(meta$disease.state)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$disease.state), col=c(1:length(levels(meta$disease.state))), pch=16, cex=0.8)

plot(mds$points, col=as.numeric(as.factor(meta$sex)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$sex), col=c(1:length(levels(meta$sex))), pch=16, cex=0.8)

plot(mds$points, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$batch), col=c(1:length(levels(meta$batch))), pch=16, cex=0.8)

plot(mds$points, col=as.numeric(as.factor(meta$tissue)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$tissue), col=c(1:length(levels(meta$tissue))), pch=16, cex=0.8)
}


#DGE ####
#function####
DGE_SEX <- function(fit, contrasts, genesOfInterest) {
  ALL_DGE <- vector(mode = "list", length = length(contrasts))
  names(ALL_DGE) <- contrasts
  DGE <- list()
  hits <- integer()
  for(i in 1:length(contrasts)) {
    topT <- topTable(fit, coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    full <- length(topT$logFC)
    topT <- topTable(fit[rownames(fit) %in% genesOfInterest$PROBEID,], coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    red <- length(topT$logFC)
    print(paste0(i, ": ", full, "/", red))
    if(red>0){
      DGE <- c(DGE, ALL_DGE[i])
      hits <- c(hits, i)
    }
  }
  
  count <- 0
  for(i in hits){
    count <- count + 1
    topT <- topTable(fit[rownames(fit) %in% genesOfInterest$PROBEID,], coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    topT[1:3,]
    names <- AnnotationDbi::select(hgu133a.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
    topT$name <- names$SYMBOL[match(rownames(topT), names$PROBEID)]
    topT <- topT[c("logFC", "AveExpr", "adj.P.Val", "name")]
    topT <- topT[order(abs(topT$logFC), decreasing = T),]
    topT <- topT[!duplicated(topT$name),]
    
    if(grepl("_F", names(DGE[count]))) {
      topT$sex <- "F"
    } else {
      topT$sex <- "M"
    }
    
    DGE[[count]] <- topT
  }
  return(DGE)
}

DGE_REG <- function(fit, contrasts, genesOfInterest) {
  ALL_DGE <- vector(mode = "list", length = length(contrasts))
  names(ALL_DGE) <- contrasts
  DGE <- list()
  hits <- integer()
  for(i in 1:length(contrasts)) {
    topT <- topTable(fit, coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    full <- length(topT$logFC)
    topT <- topTable(fit[rownames(fit) %in% genesOfInterest$PROBEID,], coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    red <- length(topT$logFC)
    print(paste0(i, ": ", full, "/", red))
    if(red>0){
      DGE <- c(DGE, ALL_DGE[i])
      hits <- c(hits, i)
    }
  }
  
  count <- 0
  for(i in hits){
    count <- count + 1
    topT <- topTable(fit[rownames(fit) %in% genesOfInterest$PROBEID,], coef = i, number = Inf)
    topT <- topT[topT$adj.P.Val < .05,]
    topT[1:3,]
    names <- AnnotationDbi::select(hgu133a.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
    topT$name <- names$SYMBOL[match(rownames(topT), names$PROBEID)]
    topT <- topT[c("logFC", "AveExpr", "adj.P.Val", "name")]
    topT <- topT[order(abs(topT$logFC), decreasing = T),]
    topT <- topT[!duplicated(topT$name),]
    
    topT$sex <- "B"
    
    DGE[[count]] <- topT
  }
  return(DGE)
}

#annot####
library(hgu133a.db)
# load("./working_data/ecs_tables.RData")
# keys <- c(cholinergic_core_genes$ensg[!cholinergic_core_genes$group %in% c("ffar")], cholinergic_tfs$tf_ensg)
keys <- readRDS("working_data/subset_keys_ensg.rds")

genesOfInterest <- select(hgu133a.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- unique(genesOfInterest)

#limma####
library(limma)
eset <- readRDS("./working_data/subset/bahn_regressed_combat.rds")
meta <- pData(eset)
plotMD(eset)

design <- model.matrix(~ 0 + disease.state + sex + hemi + tissue, data = meta)
colnames(design) <- c("BD", "CTL", "sexM", "hemiR", "tissueOFC")
contrasts <- c("BD-CTL")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennDiagram(results)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_REG(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/bahn_chol_bd.rds")


#sex effect?####
levels(meta$sex) <- c("F", "M")
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))

design <- model.matrix(~ 0 + dis_sex + tissue + hemi, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
colnames(design) <- gsub(" ", "", colnames(design))
contrasts <- c("BD_F-CTL_F", "BD_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(eset, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennDiagram(results)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/bahn_chol_bd_sex.rds")

#MIR TARGETS####

library(hgu133a.db)
mir_targets <- readRDS("./working_data/mirna_seq/mir_sequencing_top_gene_targets_weighted.rds")
keys <- mir_targets$g.ensg[mir_targets$top]
genesOfInterest <- select(hgu133a.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- unique(genesOfInterest)

#limma####
library(limma)
plotMD(eset)

design <- model.matrix(~ 0 + disease.state + sex + tissue + hemi, data = meta)
colnames(design) <- gsub("disease.state", "", colnames(design))
colnames(design) <- gsub(" ", "", colnames(design))
contrasts <- c("BD-CTL")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennDiagram(results)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_REG(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/bahn_mtar_bd.rds")


#sex effect?####
levels(meta$sex) <- c("F", "M")
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))

design <- model.matrix(~ 0 + dis_sex + tissue + hemi, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
colnames(design) <- gsub(" ", "", colnames(design))
contrasts <- c("BD_F-CTL_F", "BD_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(eset, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennDiagram(results)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_REG(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/bahn_mtar_bd_sex.rds")
