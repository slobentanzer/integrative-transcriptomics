rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#GSE53987

library(BiocParallel)
registered()
library(lme4)
library(GEOquery)
library(Biobase)
library(BiocParallel)
library(affy)
library(sva)
library(hugene10stv1cdf)
library(limma)
library(ggplot2)


if(!file.exists("./working_data/subset/lanz_raw_eset.rds")) {
  #CEL READ####
  affy.folder <- "./raw_data/datasets/Lanz_GSE53987/"
  affy.files <- paste0(affy.folder, list.files(affy.folder))
  affy.files <- affy.files[grep("CEL", affy.files)]
  
  library(affy)
  raw_read <- ReadAffy(filenames = affy.files)
  eset <- rma(raw_read)
  eset_exprs <- exprs(eset)
  names <- colnames(eset_exprs)
  names <- unlist(lapply(strsplit(names, "_", fixed = T), "[", 1))
  colnames(eset_exprs) <- names
  limma::plotMA(eset_exprs)
  
  #META####
  if (!file.exists("./05_subset_analysis/02_lanz/gpl570.rds")) {
    gpl570 <- getGEO("GSE53987")
    gpl570 <- gpl570[[1]]
    save(gpl570, file = "./05_subset_analysis/02_lanz/gpl570.rds")
  } else {
    load("./05_subset_analysis/02_lanz/gpl570.rds")
  }
  # exprs(gpl570) <- eset_exprs
  
  #format gpl570####
  pheno <- pData(gpl570)
  head(pheno)
  colnames(eset_exprs) == rownames(pheno)
  disease.state <- pheno$characteristics_ch1.7
  disease.state <-
    strsplit(as.character(disease.state), ": ", fixed = T)
  disease.state <- unlist(lapply(disease.state, "[", 2))
  disease.state[grep("bipol", disease.state)] <- "BD"
  disease.state[grep("schiz", disease.state)] <- "SCZ"
  disease.state[grep("major", disease.state)] <- "MDD"
  disease.state[grep("contr", disease.state)] <- "CTL"
  disease.state <- factor(disease.state)
  
  region <- pheno$characteristics_ch1.6
  region <- strsplit(as.character(region), ": ", fixed = T)
  region <- unlist(lapply(region, "[", 2))
  region[grep("hippo", region)] <- "HIP"
  region[grep("front", region)] <- "PFC"
  region[grep("stria", region)] <- "STR"
  region <- factor(region)
  
  reg_dis <- factor(paste(region, disease.state, sep = "_"))
  
  #confounders
  age <- pheno$characteristics_ch1
  age <- strsplit(as.character(age), ": ", fixed = T)
  age <- as.numeric(unlist(lapply(age, "[", 2)))
  sex <- pheno$characteristics_ch1.1
  sex <- strsplit(as.character(sex), ": ", fixed = T)
  sex <- factor(unlist(lapply(sex, "[", 2)))
  race <- pheno$characteristics_ch1.2
  race <- strsplit(as.character(race), ": ", fixed = T)
  race <- factor(unlist(lapply(race, "[", 2)))
  pmi <- pheno$characteristics_ch1.3
  pmi <- strsplit(as.character(pmi), ": ", fixed = T)
  pmi <- as.numeric(unlist(lapply(pmi, "[", 2)))
  ph <- pheno$characteristics_ch1.4
  ph <- strsplit(as.character(ph), ": ", fixed = T)
  ph <- as.numeric(unlist(lapply(ph, "[", 2)))
  rin <- pheno$characteristics_ch1.5
  rin <- strsplit(as.character(rin), ": ", fixed = T)
  rin <- as.numeric(unlist(lapply(rin, "[", 2)))
  prot <- protocolData(raw_read)
  batch <- factor(substring(prot$ScanDate, 1, nchar(prot$ScanDate) - 9))
  
  reg_dis_sex <- factor(paste(reg_dis, sex, sep = "_"))
  
  RNAdeg <- AffyRNAdeg(raw_read)
  rna.deg <- RNAdeg$slope
  
  meta <-
    data.frame(
      disease.state,
      region,
      reg_dis,
      reg_dis_sex,
      age,
      sex,
      race,
      pmi,
      ph,
      rin,
      batch,
      rna.deg
    )
  rownames(meta) <- colnames(eset_exprs)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  eset
  
  saveRDS(eset, "./working_data/subset/lanz_raw_eset.rds")
}
#QC, correlation####

if(!file.exists("./working_data/subset/lanz_eset_combat_regressed.rds")) {
  eset <- readRDS("./working_data/subset/lanz_raw_eset.rds")
  meta <- pData(eset)
  
  #batches
  table(meta$batch)
  idx <- which(meta$batch == "2012-07-31T")
  eset <- eset[,-idx]
  meta <- pData(eset)
  
  
  #groups
  plot(meta$disease.state)
  
  #colinearities
  for(i in 2:length(meta)) {
    cor <- anova(lm(as.numeric(meta[,i]) ~ meta$disease.state))
    plot(meta[,i] ~ meta$disease.state, main = paste(colnames(meta)[i], " =
                                                     ", cor$`Pr(>F)`[1]))
    if(cor$`Pr(>F)`[1] < .05)
      print(colnames(meta)[i])
  }
  plot(meta$race)
  
  hist(meta$ph)
  
  #CORRECTION####
  #batch
  meta$batch <- factor(meta$batch)
  library(sva)
  model <- model.matrix(~disease.state, data = meta)
  eset_exprs <- ComBat(exprs(eset), meta$batch, model)
  
  #regression
  X <- model.matrix(~disease.state+region+sex+age+pmi+ph+race+rin+rna.deg, data=meta)
  Y <- eset_exprs
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  reg_corr <- (as.matrix(X[,9:13]) %*% (as.matrix(beta[9:13,])))
  eset_exprs <- eset_exprs - t(reg_corr)
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  saveRDS(eset, "./working_data/subset/lanz_eset_combat_regressed.rds")
}
if(F){
eset <- readRDS("./working_data/subset/lanz_eset_combat_regressed.rds")
meta <- pData(eset)

#MDS
mds <- cmdscale(dist(t(exprs(eset))), eig = T)
plot(mds$points)
pc1 <- mds$eig[1]^2 / sum(mds$eig^2)
pc2 <- mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(meta$disease.state)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$disease.state), col=c(1:length(levels(meta$disease.state))), pch=16, cex=0.8)

plot(mds$points, col=as.numeric(as.factor(meta$region)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$region), col=c(1:length(levels(meta$region))), pch=16, cex=0.8)

plot(mds$points, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot", asp=1, 
     xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
legend("bottomright", levels(meta$batch), col=c(1:length(levels(meta$batch))), pch=16, cex=0.8)
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


#ANNOTATION####
# [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]

library(hgu133plus2.db)

# load("./working_data/ecs_tables.RData")
# keys <- c(cholinergic_core_genes$ensg, cholinergic_tfs$tf_ensg)
keys <- readRDS("working_data/subset_keys_ensg.rds")

genesOfInterest <- AnnotationDbi::select(hgu133plus2.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[order(genesOfInterest$PROBEID),]
genesOfInterest <- unique(genesOfInterest)
genesOfInterest <- genesOfInterest[-which(duplicated(genesOfInterest$PROBEID)),]

#limma####
eset <- readRDS("./working_data/subset/lanz_eset_combat_regressed.rds")
meta <- pData(eset)

library(limma)

design <- model.matrix(~ 0 + reg_dis + age + sex + pmi + rin + rna.deg, data = meta)
colnames(design) <- gsub("reg_dis", "", colnames(design))
contrasts <- c("HIP_SCZ-HIP_CTL", "HIP_BD-HIP_CTL",
               "PFC_SCZ-PFC_CTL", "PFC_BD-PFC_CTL",
               "STR_SCZ-STR_CTL", "STR_BD-STR_CTL")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
# vennCounts(results)

DGE <- DGE_REG(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/lanz_chol_bd_scz.rds")

#SEX DIFF####
design <- model.matrix(~ 0 + reg_dis_sex + age + pmi + rin + rna.deg, data = meta)
colnames(design) <- gsub("reg_dis_sex", "", colnames(design))
contrasts <- c("HIP_SCZ_F-HIP_CTL_F", "HIP_BD_F-HIP_CTL_F",
"PFC_SCZ_F-PFC_CTL_F", "PFC_BD_F-PFC_CTL_F",
"STR_SCZ_F-STR_CTL_F", "STR_BD_F-STR_CTL_F", 
"HIP_SCZ_M-HIP_CTL_M", "HIP_BD_M-HIP_CTL_M",
"PFC_SCZ_M-PFC_CTL_M", "PFC_BD_M-PFC_CTL_M", 
"STR_SCZ_M-STR_CTL_M", "STR_BD_M-STR_CTL_M")

contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
# vennCounts(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/lanz_chol_sex_bd_scz.rds")

#MIR TARGETS####
mir_targets <- readRDS("./working_data/mirna_seq/mir_sequencing_top_gene_targets_weighted.rds")
keys <- c(mir_targets$g.ensg[mir_targets$top])

genesOfInterest <- AnnotationDbi::select(hgu133plus2.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[order(genesOfInterest$PROBEID),]
genesOfInterest <- unique(genesOfInterest)

#limma####
library(limma)
plotMD(eset)

design <- model.matrix(~ 0 + reg_dis + age + sex + pmi + rin + rna.deg, data = meta)
colnames(design) <- gsub("reg_dis", "", colnames(design))
contrasts <- c("HIP_SCZ-HIP_CTL", "HIP_BD-HIP_CTL",
               "PFC_SCZ-PFC_CTL", "PFC_BD-PFC_CTL",
               "STR_SCZ-STR_CTL", "STR_BD-STR_CTL")
contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

fit <- lmFit(eset, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
# vennCounts(results)

DGE <- DGE_REG(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/lanz_mtar_bd_scz.rds")

#SEX DIFF####
design <- model.matrix(~ 0 + reg_dis_sex + age + pmi + rin + rna.deg, data = meta)
colnames(design) <- gsub("reg_dis_sex", "", colnames(design))
contrasts <- c("HIP_SCZ_F-HIP_CTL_F", "HIP_BD_F-HIP_CTL_F", 
               "PFC_SCZ_F-PFC_CTL_F", "PFC_BD_F-PFC_CTL_F", 
               "STR_SCZ_F-STR_CTL_F", "STR_BD_F-STR_CTL_F", 
               "HIP_SCZ_M-HIP_CTL_M", "HIP_BD_M-HIP_CTL_M", 
               "PFC_SCZ_M-PFC_CTL_M", "PFC_BD_M-PFC_CTL_M", 
               "STR_SCZ_M-STR_CTL_M", "STR_BD_M-STR_CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts, 
                                 levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
# vennCounts(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

#save####
saveRDS(DGE, file = "./working_data/subset/dge/lanz_mtar_bd_scz_sex.rds")
