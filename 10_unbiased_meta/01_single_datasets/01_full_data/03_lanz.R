#LOAD LANZ DATASET, NORMALIZE, BATCH CORRECT, REGRESS COVARIATES####
#GSE53987####

rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#01 read data####
#meta####
meta = read.csv("./raw_data/datasets/Lanz_GSE53987/Lanz_GSE53987_datMeta.csv")
rownames(meta) = meta$Chip
meta$Sex = as.factor(meta$Sex)
meta$Race = as.factor(meta$Race)
meta$Region[meta$BrainRegion == "hippocampus"] = "HC"
meta$Region[meta$BrainRegion == "Pre-frontal cortex (BA46)"] = "PFC"
meta$Region[meta$BrainRegion == "Associative striatum"] = "STR"
meta$Region = as.factor(meta$Region)
meta$Group[meta$Disorder == "bipolar disorder"] = "BD"
meta$Group[meta$Disorder == "control"] = "CTL"
meta$Group[meta$Disorder == "major depressive disorder"] = "MDD"
meta$Group[meta$Disorder == "schizophrenia"] = "SCZ"
meta$Group = factor(meta$Group, levels = c("CTL", "BD", "SCZ", "MDD"))
meta$Study = "SCZ.BD.MDD.Lanz"
meta$Subject = paste(
  meta$BrainBank,
  meta$Group,
  meta$Sex,
  meta$Age,
  meta$PMI,
  meta$Race,
  sep = "."
)

#exprs####
pfc_samples = list.files("./raw_data/datasets/Lanz_GSE53987/")
pfc_samples = pfc_samples[grep("B46", pfc_samples)]
data.affy = ReadAffy(celfile.path = "./raw_data/datasets/Lanz_GSE53987/", filenames =
                       pfc_samples)
exprs = exprs(data.affy)

samples = substr(colnames(exprs), 1, 10)
meta = meta[samples, ]
sd = protocolData(data.affy)$ScanDate
sd = substring(sd, 1, nchar(sd) - 9)
meta$Batch = sd

#RNAdeg
RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
meta$RNAdeg = RNAdeg$slope

#rma
exprs = rma(data.affy, background = T, normalize = T)
exprs = exprs(exprs)
colnames(exprs) = substr(colnames(exprs), 1, 10)

#annotation (gandal)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host =
                    "feb2014.archive.ensembl.org")
identifier <- "affy_hg_u133_plus_2"
getinfo <-
  c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "entrezgene",
    "external_gene_id",
    "chromosome_name",
    "start_position",
    "end_position"
  )
geneDat <-
  getBM(
    attributes = getinfo,
    filters = identifier,
    values = rownames(exprs),
    mart = ensembl
  )

idx = match(rownames(exprs), geneDat$affy_hg_u133_plus_2)
annot = cbind(rownames(exprs), geneDat[idx, ])

save(file = "./working_data/01_unbiased/lanz_compiled.RData", exprs, meta, annot)


#02 normalization and regression####
load("./working_data/01_unbiased/lanz_compiled.RData")

#batch control
table(meta$Batch)

#remove outliers (gandal)
sdout <- 2
normadj <- (0.5+0.5*cor(exprs, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
print(colnames(exprs)[outliers])
print(table(outliers))
plot(z.ku, col = as.numeric(meta$Group), pch=19, main="Outlier detection", ylab="Euclidean distance (z score)")
abline(h=-2, lty=2)
exprs = exprs[,!outliers]
meta = meta[!outliers,]

#sub quality (gandal)
to_remove = (meta$Group=="SCZ") & (meta$RNAdeg > 5)
exprs = exprs[,!to_remove]
meta = meta[!to_remove,]

#combat
mod = model.matrix(~meta$Group)
batch = as.factor(meta$Batch)
exprs.combat = ComBat(exprs, batch, mod)
exprs = exprs.combat

#collapse rows
realGenes = !is.na(annot$ensembl_gene_id)  #
exprs = exprs[realGenes,]
annot = annot[realGenes,]

CR = collapseRows(exprs, rowGroup = annot$ensembl_gene_id, rowID = annot$affy_hg_u133_plus_2) 
exprs = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], annot$affy_hg_u133_plus_2)
annot = annot[idx,]
rownames(annot) = annot$ensembl_gene_id


#covariate regression
X = model.matrix(~Group+Sex+Age+PMI+pH+Race+RIN+RNAdeg, data=meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,5:11]) %*% (as.matrix(beta[5:11,])))  
exprs = exprs - t(to_regress)

#randomize ctl
idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/01_unbiased/lanz_normalized.RData", exprs, meta, annot)

  
  