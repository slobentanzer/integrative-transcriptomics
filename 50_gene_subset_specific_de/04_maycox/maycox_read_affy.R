rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#GSE17612

library(BiocParallel)
registered()
library(lme4)
library(GEOquery)
library(Biobase)
library(BiocParallel)
library(affy)
library(hugene10stv1cdf)
library(limma)


if(!file.exists("./working_data/subset/maycox_raw_eset.rds")) {
  
  #READ AFFY####
  library(affy)
  affy.folder <- "./raw_data/datasets/Maycox_GSE17612/"
  affy.files <- paste0(affy.folder, list.files(affy.folder))
  affy.files <- affy.files[grep("CEL", affy.files)]
  raw_read <- ReadAffy(filenames = affy.files)
  eset <- rma(raw_read)
  
  meta <-
    read.csv("./05_subset_analysis/04_maycox/Maycox_GSE17612_datMeta.csv")
  rownames(meta) = meta$Chip
  colnames(meta)[6:9] <- tolower(colnames(meta)[6:9])
  meta$disease.state = factor(meta$Group, levels = c("CTL", "SCZ"))
  meta$sex = as.factor(meta$sex)
  head(meta)
  
  eset_exprs <- exprs(eset)
  names <- colnames(eset_exprs)
  names <- unlist(lapply(strsplit(names, ".", fixed = T), "[", 1))
  colnames(eset_exprs) <- names
  plotMA(eset_exprs)
  
  idx <- match(colnames(eset_exprs), rownames(meta))
  meta = meta[idx, ]
  
  RNAdeg = AffyRNAdeg(raw_read)
  plotAffyRNAdeg(RNAdeg)
  meta$rna.deg = RNAdeg$slope
  
  batch <- factor(substring(protocolData(raw_read)$ScanDate, 1, 8))
  batch <- data.frame(names, batch)
  
  idx <- match(rownames(meta), batch$names)
  meta$batch <- batch$batch[idx]
  meta$batch <- factor(meta$batch)
  
  mds <- cmdscale(dist(t(eset_exprs)))
  plot(mds,
       col = as.numeric(as.factor(meta$disease.state)),
       pch = 16,
       main = "MDS Plot")
  legend(
    "topright",
    levels(meta$disease.state),
    col = c(1:2),
    pch = 16,
    cex = 0.8
  )
  
  plot(mds,
       col = as.numeric(as.factor(meta$batch)),
       pch = 16,
       main = "MDS Plot")
  legend(
    "topright",
    levels(meta$batch),
    col = seq(length(levels(meta$batch))),
    pch = 16,
    cex = 0.8
  )
  
  plot(mds,
       col = as.numeric(as.factor(meta$sex)),
       pch = 16,
       main = "MDS Plot")
  legend(
    "topright",
    levels(meta$sex),
    col = seq(length(levels(meta$sex))),
    pch = 16,
    cex = 0.8
  )
  
  rownames(meta) == colnames(eset_exprs)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  eset
  
  saveRDS(eset, "./working_data/subset/maycox_raw_eset.rds")
}
#QC, correlation####
if(!file.exists("./working_data/subset/maycox_regressed_ph_rna_combat.rds")) {
  eset <- readRDS("./working_data/subset/maycox_raw_eset.rds")
  meta <- pData(eset)[,c("disease.state", "sex", "age", "pmi", "ph", "rna.deg", "batch")]
  
  #batches
  table(meta$batch)
  
  #groups
  plot(meta$disease.state)
  head(meta)
  
  #colinearities
  for(i in 2:length(meta)) {
    cor <- anova(lm(as.numeric(meta[,i]) ~ meta$disease.state))
    plot(meta[,i] ~ meta$disease.state, main = paste(colnames(meta)[i], " =
                                                     ", cor$`Pr(>F)`[1]))
    if(cor$`Pr(>F)`[1] < .05)
      print(colnames(meta)[i])
  }
  hist(meta$ph)
  
  #CORRECTION####
  library(sva)
  meta <- pData(eset)
  meta$batch <- factor(meta$batch)
  meta$disease.state <- factor(meta$disease.state)
  plot(meta$disease.state ~ meta$batch)
  model = model.matrix(~disease.state, data=meta)
  eset_exprs = ComBat(exprs(eset), batch=meta$batch, mod=model)
  
  par(mfrow = c(1,2))
  mds <- cmdscale(dist(t(exprs(eset))))
  plot(mds, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot")
  legend("topleft", levels(meta$batch), col=c(1:3), pch=16, cex=0.8)
  
  mds.combat <- cmdscale(dist(t(eset_exprs)))
  plot(mds.combat, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot ComBat")
  legend("topleft", levels(meta$batch), col=c(1:3), pch=16, cex=0.8)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  #regression, ph, rna.deg
  
  eset_exprs <- exprs(eset)
  X <- model.matrix(~disease.state+sex+age+ph+pmi+rna.deg, data=meta)
  Y <- eset_exprs
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  reg_corr <- (as.matrix(X[,c("ph", "rna.deg")]) %*% (as.matrix(beta[c("ph", "rna.deg"),])))
  eset_exprs <- eset_exprs - t(reg_corr)
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  saveRDS(eset, "./working_data/subset/maycox_regressed_ph_rna_combat.rds")
}
1
if(F){
eset <- readRDS("./working_data/subset/maycox_regressed_ph_rna_combat.rds")
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

plot(mds$points, col=as.numeric(as.factor(meta$ph)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(as.factor(meta$ph)), col=c(1:length(levels(as.factor(meta$ph)))), pch=16, cex=0.8)
}

#DGE####
#function####
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
    mean(genesOfInterest$PROBEID[which(genesOfInterest$ENSEMBL %in% keys)] %in% rownames(topT))
    topT <- topT[topT$adj.P.Val < .05,]
    topT[1:3,]
    names <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
    topT$name <- names$SYMBOL[match(rownames(topT), names$PROBEID)]
    topT <- topT[c("logFC", "AveExpr", "adj.P.Val", "name")]
    topT <- topT[order(abs(topT$logFC), decreasing = T),]
    topT <- topT[!duplicated(topT$name),]
    
    topT$sex <- "B"
    
    DGE[[count]] <- topT
  }
  return(DGE)
}

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
    mean(genesOfInterest$PROBEID[which(genesOfInterest$ENSEMBL %in% keys)] %in% rownames(topT))
    topT <- topT[topT$adj.P.Val < .05,]
    topT[1:3,]
    names <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
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



#U133 Plus####
library(hgu133plus2.db)
# load("./working_data/ecs_tables.RData")
# keys <- c(cholinergic_core_genes$ensg[!cholinergic_core_genes$group %in% c("ffar")], cholinergic_tfs$tf_ensg)
keys <- readRDS("working_data/subset_keys_ensg.rds")

genesOfInterest <- select(hgu133plus2.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- unique(genesOfInterest)


#limma####
library(limma)
eset <- readRDS("./working_data/subset/maycox_regressed_ph_rna_combat.rds")
meta <- pData(eset)
plotMD(eset)

design <- model.matrix(~ 0 + disease.state + sex + pmi + age, data = meta)
colnames(design) <- gsub("disease.state", "", colnames(design))
contrasts <- c("SCZ-CTL")
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
saveRDS(DGE, file = "./working_data/subset/dge/maycox_chol_scz.rds")


#SEX DIFF####
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))
design <- model.matrix(~ 0 + dis_sex + pmi + age, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
contrasts <- c("SCZ_F-CTL_F", 
               "SCZ_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts,
                                 levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

saveRDS(DGE, file = "./working_data/subset/dge/maycox_chol_scz_sex.rds")


#MIR TARGETS####
library(hgu133plus2.db)
mir_targets <- readRDS("./working_data/mirna_seq/mir_sequencing_top_gene_targets_weighted.rds")
keys <- mir_targets$g.ensg[mir_targets$top]
genesOfInterest <- select(hgu133plus2.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- unique(genesOfInterest)


#limma####
library(limma)
plotMD(eset)

meta <- pData(eset)

design <- model.matrix(~ 0 + disease.state + sex + pmi + age, data = meta)
colnames(design) <- gsub("disease.state", "", colnames(design))
contrasts <- c("SCZ-CTL")
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
saveRDS(DGE, file = "./working_data/subset/dge/maycox_mtar_scz.rds")


#SEX DIFF####
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))
design <- model.matrix(~ 0 + dis_sex + pmi + age, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
contrasts <- c("SCZ_F-CTL_F", 
               "SCZ_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts,
                                 levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit[rownames(fit) %in% genesOfInterest$PROBEID,])
vennDiagram(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

saveRDS(DGE, file = "./working_data/subset/dge/maycox_mtar_scz_sex.rds")

