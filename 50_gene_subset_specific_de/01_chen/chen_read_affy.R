rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#GSE35978

library(BiocParallel)
registered()
library(lme4)
library(GEOquery)
library(Biobase)
library(BiocParallel)
library(affy)
library(sva)
# library(hugene10stv1cdf)
library(limma)


if(!file.exists("working_data/subset/chen_raw_eset.rds")){
  #CEL READ####
  affy.folder <- "./raw_data/datasets/Chen_GSE35978/"
  affy.files <- paste0(affy.folder, list.files(affy.folder))
  affy.files <- affy.files[grep("CEL", affy.files)]
  
  raw_read <- ReadAffy(filenames = affy.files)
  eset <- rma(raw_read)
  eset_exprs <- exprs(eset)
  names <- colnames(eset_exprs)
  names <- unlist(lapply(strsplit(names, ".", fixed = T), "[", 1))
  colnames(eset_exprs) <- names
  limma::plotMA(eset_exprs)
  
  batch <- factor(substring(protocolData(raw_read)$ScanDate, 1, 10))
  batch <- data.frame(names, batch)
  
  #META####
  if (!file.exists("./05_subset_analysis/01_chen/gpl6244.rds"))
  {
    gpl6244 <- getGEO("GSE35978")
    gpl6244 <- gpl6244[[1]]
    saveRDS(gpl6244, file = "./05_subset_analysis/01_chen/gpl6244.rds")
  } else
    gpl6244 <- readRDS(file = "./05_subset_analysis/01_chen/gpl6244.rds")
  
  pheno <- pData(gpl6244)
  rownames(pheno) == colnames(eset_exprs)
  
  head(pheno)
  tissue <- pheno$characteristics_ch1
  tissue <- strsplit(as.character(tissue), ": ", fixed = T)
  tissue <- unlist(lapply(tissue, "[", 2))
  tissue[grep("pari", tissue)] <- "CTXp"
  tissue[grep("cereb", tissue)] <- "CER"
  tissue <- factor(tissue)
  
  disease.state <- pheno$characteristics_ch1.1
  disease.state <- strsplit(as.character(disease.state), ": ", fixed = T)
  disease.state <- unlist(lapply(disease.state, "[", 2))
  disease.state[grep("schiz", disease.state)] <- "SCZ"
  disease.state[grep("bipo", disease.state)] <- "BD"
  disease.state[grep("depr", disease.state)] <- "MDD"
  disease.state[grep("unaff", disease.state)] <- "CTL"
  disease.state <- factor(disease.state)
  
  tis_dis <- factor(paste(tissue, disease.state, sep = "_"))
  
  age <- pheno$characteristics_ch1.2
  age <- strsplit(as.character(age), ": ", fixed = T)
  age <- as.numeric(unlist(lapply(age, "[", 2)))
  ph <- pheno$characteristics_ch1.5
  ph <- strsplit(as.character(ph), ": ", fixed = T)
  ph <- as.numeric(unlist(lapply(ph, "[", 2)))
  sex <- pheno$characteristics_ch1.4
  sex <- strsplit(as.character(sex), ": ", fixed = T)
  sex <- factor(unlist(lapply(sex, "[", 2)))
  pmi <- pheno$characteristics_ch1.3
  pmi <- strsplit(as.character(pmi), ": ", fixed = T)
  pmi <- as.numeric(unlist(lapply(pmi, "[", 2)))
  
  RNAdeg <- AffyRNAdeg(raw_read)
  rna.deg <- RNAdeg$slope
  
  meta <- data.frame(disease.state, tissue, tis_dis, age, ph, sex, pmi, rna.deg)
  rownames(meta) <- rownames(pheno)
  
  idx <- match(rownames(meta), batch$names)
  meta$batch <- batch$batch[idx]
  meta$batch <- factor(meta$batch)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  saveRDS(eset, "./working_data/subset/chen_raw_eset.rds")
}

if(!file.exists("./working_data/subset/chen_eset_ctx_normalized.rds")){
  eset <- readRDS("./working_data/subset/chen_raw_eset.rds")
  
  #remove na and cerebellum
  idx <- which(meta$disease.state == "NA" | meta$tissue == "CER")
  meta <- meta[-idx,]
  eset_exprs <- eset_exprs[,-idx]
  meta$disease.state <- factor(meta$disease.state)
  meta$tissue <- factor(meta$tissue)
  meta$batch <- factor(meta$batch)
  
  rownames(meta) == colnames(eset_exprs)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  eset
  
  saveRDS(eset, "./working_data/subset/chen_eset_ctx_normalized.rds")
}
#QC, correlation####
if(!file.exists("./working_data/subset/chen_eset_regressed_ph_pmi_combat.rds")){
  eset <- readRDS("./working_data/subset/chen_eset_ctx_normalized.rds")
  meta <- pData(eset)
  
  #batches
  table(meta$batch)
  idx <- which(meta$batch == "2010-04-14")
  meta[idx,]
  eset <- eset[,-idx]
  meta <- pData(eset)
  meta$batch <- factor(meta$batch)
  meta$sex <- factor(meta$sex)
  
  #groups
  plot(meta$disease.state)
  head(meta)
  
  #colinearities
  for(i in 4:length(meta)) {
    cor <- anova(lm(as.numeric(meta[,i]) ~ meta$disease.state))
    plot(meta[,i] ~ meta$disease.state, main = paste(colnames(meta)[i], " =
                                                     ", cor$`Pr(>F)`[1]))
    if(cor$`Pr(>F)`[1] < .05)
      print(colnames(meta)[i])
  }
  hist(meta$ph)
  hist(meta$pmi)
  plot(meta$batch)
  
  #CORRECTION####
  #batch
  plot(meta$disease.state ~ meta$batch)
  model = model.matrix(~disease.state, data=meta)
  eset_exprs = ComBat(exprs(eset), batch=meta$batch, mod=model)
  
  mds <- cmdscale(dist(t(exprs(eset))))
  plot(mds, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot")
  legend("topright", levels(meta$batch), col=c(1:3), pch=16, cex=0.8)
  
  mds.combat <- cmdscale(dist(t(eset_exprs)))
  plot(mds.combat, col=as.numeric(as.factor(meta$batch)), pch=16, main="MDS Plot ComBat")
  legend("topright", levels(meta$batch), col=c(1:3), pch=16, cex=0.8)
  
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  #regression, ph, pmi
  eset_exprs <- exprs(eset)
  X <- model.matrix(~disease.state+age+sex+pmi+ph+rna.deg, data=meta)
  Y <- eset_exprs
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  
  reg_corr <- (as.matrix(X[,7:9]) %*% (as.matrix(beta[7:9,])))
  eset_exprs <- eset_exprs - t(reg_corr)
  eset <- ExpressionSet(eset_exprs)
  pData(eset) <- meta
  
  saveRDS(eset, "./working_data/subset/chen_eset_regressed_ph_pmi_combat.rds")
}
1

if(F){
  
  eset <-
    readRDS("./working_data/subset/chen_eset_regressed_ph_pmi_combat.rds")
  meta <- pData(eset)
  
  #MDS
  mds <- cmdscale(dist(t(exprs(eset))), eig = T)
  plot(mds$points)
  pc1 <- mds$eig[1] ^ 2 / sum(mds$eig ^ 2)
  pc2 <- mds$eig[2] ^ 2 / sum(mds$eig ^ 2)
  plot(
    mds$points,
    col = as.numeric(as.factor(meta$disease.state)),
    pch = 16,
    main = "MDS Plot",
    asp = 1,
    xlab = paste("PC1 (", signif(100 * pc1, 3), "%)", sep = ""),
    ylab = paste("PC1 (", signif(100 * pc2, 3), "%)", sep = "")
  )
  legend(
    "bottomright",
    levels(meta$disease.state),
    col = c(1:length(levels(meta$disease.state))),
    pch = 16,
    cex = 0.8
  )
  
  plot(
    mds$points,
    col = as.numeric(as.factor(meta$sex)),
    pch = 16,
    main = "MDS Plot",
    asp = 1,
    xlab = paste("PC1 (", signif(100 * pc1, 3), "%)", sep = ""),
    ylab = paste("PC1 (", signif(100 * pc2, 3), "%)", sep = "")
  )
  legend(
    "bottomright",
    levels(meta$sex),
    col = c(1:length(levels(meta$sex))),
    pch = 16,
    cex = 0.8
  )
  
  plot(
    mds$points,
    col = as.numeric(as.factor(meta$batch)),
    pch = 16,
    main = "MDS Plot",
    asp = 1,
    xlab = paste("PC1 (", signif(100 * pc1, 3), "%)", sep = ""),
    ylab = paste("PC1 (", signif(100 * pc2, 3), "%)", sep = "")
  )
  legend(
    "bottomright",
    levels(meta$batch),
    col = c(1:length(levels(meta$batch))),
    pch = 16,
    cex = 0.8
  )
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
    names <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
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
    names <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = as.character(rownames(topT)), columns = c("SYMBOL"), keytype = "PROBEID")
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

#annot####
library(hugene10sttranscriptcluster.db)

# load("./working_data/ecs_tables.RData")
# keys <- c(cholinergic_core_genes$ensg[!cholinergic_core_genes$group %in% c("ffar")], cholinergic_tfs$tf_ensg)
keys <- readRDS("working_data/subset_keys_ensg.rds")

genesOfInterest <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[!is.na(genesOfInterest$SYMBOL),]


#limma####
eset <- readRDS("./working_data/subset/chen_eset_regressed_ph_pmi_combat.rds")
meta <- pData(eset)

library(limma)
plotMD(exprs(eset))

design <- model.matrix(~ 0 + disease.state + age + sex, data = meta)
colnames(design) <- gsub("disease.state", "", colnames(design))
colnames(design) <- gsub("BP", "BD", colnames(design))
contrasts <- c("SCZ-CTL", "BD-CTL")
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
saveRDS(DGE, file = "./working_data/subset/dge/chen_chol_scz.rds")


#SEX DIFF####
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))
design <- model.matrix(~ 0 + dis_sex, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
colnames(design) <- gsub("BP", "BD", colnames(design))
contrasts <- c("SCZ_F-CTL_F", "BD_F-CTL_F",
               "SCZ_M-CTL_M", "BD_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts,
                                 levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennCounts(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

saveRDS(DGE, file = "./working_data/subset/dge/chen_chol_bd_scz_sex.rds")

#MIR TARGETS####
#annot####
library(hugene10sttranscriptcluster.db)


mir_targets <- readRDS("./working_data/mirna_seq/mir_sequencing_top_gene_targets_weighted.rds")
keys <- mir_targets$g.ensg[mir_targets$top]

genesOfInterest <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys = keys, columns = c("PROBEID", "SYMBOL"), keytype = "ENSEMBL")
genesOfInterest[is.na(genesOfInterest$SYMBOL),]
genesOfInterest <- genesOfInterest[-which(is.na(genesOfInterest$SYMBOL)),]


#limma####
library(limma)
plotMD(exprs(eset))
meta <- pData(eset)

design <- model.matrix(~ 0 + disease.state + age + sex + rna.deg, data = meta)
colnames(design) <- gsub("disease.state", "", colnames(design))
colnames(design) <- gsub("BP", "BD", colnames(design))
contrasts <- c("SCZ-CTL", "BD-CTL")
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
saveRDS(DGE, file = "./working_data/subset/dge/chen_mtar_scz.rds")

#SEX DIFF####
meta$dis_sex <- factor(paste(meta$disease.state, meta$sex, sep = "_"))
design <- model.matrix(~ 0 + dis_sex, data = meta)
colnames(design) <- gsub("dis_sex", "", colnames(design))
colnames(design) <- gsub("BP", "BD", colnames(design))
contrasts <- c("SCZ_F-CTL_F", "BD_F-CTL_F",
               "SCZ_M-CTL_M", "BD_M-CTL_M")
contrast.matrix <- makeContrasts(contrasts = contrasts,
                                 levels = design)

fit <- lmFit(exprs(eset), design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

results <- decideTests(fit)
vennCounts(results)

DGE <- DGE_SEX(fit, contrasts, genesOfInterest)

saveRDS(DGE, file = "./working_data/subset/dge/chen_mtar_bd_scz_sex.rds")
