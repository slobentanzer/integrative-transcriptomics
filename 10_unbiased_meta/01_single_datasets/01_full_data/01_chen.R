#LOAD CHEN DATASET, NORMALIZE, BATCH CORRECT, REGRESS COVARIATES####
#GSE35978####

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
meta = read.csv("./raw_data/datasets/Chen_GSE35978/GSE35978_datMeta.csv")
rownames(meta) = meta$Chip
meta$Region = as.factor(meta$Region)
meta$Group[meta$Group == "Bipolar"] = "BD"
meta$Group[meta$Group == "Control"] = "CTL"
meta$Group[meta$Group == "Schizophrenia"] = "SCZ"
meta$Group[meta$Group == "Depression"] = "MDD"

meta$Group = factor(meta$Group, levels = c("CTL", "MDD", "BD", "SCZ"))
meta$Age = as.numeric(meta$Age)
meta$PMI = as.numeric(meta$PMI)
meta$pH = as.numeric(meta$pH)
meta$Sex = as.factor(meta$Sex)
meta$Study = "SCZ.BD.Chen"

#remove NA
to_keep = !is.na(meta$Group)
meta = meta[to_keep, ]

#remove other brain regions
ctx_only = rownames(meta)[meta$Region == "PCTX"]
all_samples = list.files("./raw_data/datasets/Chen_GSE35978/")
ctx_sample_files = all_samples[pmatch(ctx_only, all_samples)]

#exprs####
data.affy = ReadAffy(celfile.path = "./raw_data/datasets/Chen_GSE35978/", filenames =
                       ctx_sample_files)
exprs = affy::rma(
  data.affy,
  verbose = T,
  normalize = T,
  background = T
)
exprs = exprs(exprs)
colnames(exprs) = substring(colnames(exprs), 1, 9)

#RNAdeg
RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
RNAdeg$sample.names = substring(RNAdeg$sample.names, 1, 9)
idx =  match(rownames(meta), RNAdeg$sample.names)
meta$RNAdeg = RNAdeg$slope[idx]

#batch
batch = as.factor(substring(protocolData(data.affy)$ScanDate, 1, 10))
idx = match(rownames(meta), colnames(exprs))
meta$Batch = batch[idx]

#meta match
idx = match(colnames(exprs), rownames(meta))
meta = meta[idx, ]

#annotation (gandal)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host =
                    "feb2014.archive.ensembl.org")
a = listAttributes(ensembl)
identifier <- "affy_hugene_1_0_st_v1"
getinfo <-
  c(
    "affy_hugene_1_0_st_v1",
    "ensembl_gene_id",
    "hgnc_symbol",
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
idx = match(rownames(exprs), geneDat$affy_hugene_1_0_st_v1)
annot = cbind(rownames(exprs), geneDat[idx, ])

#save
save(file = "./working_data/01_unbiased/chen_compiled.RData", exprs, meta, annot)

#02 normalization and regression####
load("./working_data/01_unbiased/chen_compiled.RData")

#batch control
table(meta$Batch)
to_keep = meta$Batch != "2010-04-14"
exprs = exprs[,to_keep]
meta = meta[to_keep,]

#remove MDD
to_keep = meta$Group != "MDD"
exprs = exprs[,to_keep]
meta = meta[to_keep,]

#remove sub quality (gandal)
to_keep = meta$PMI < 70
to_keep = to_keep & ((meta$pH < 6.9) | !(meta$Group=="CTL")) & (meta$pH > 5.8)
table(to_keep)
exprs = exprs[,to_keep]
meta = meta[to_keep,]
meta$Group = factor(meta$Group)

#remove outliers (gandal)
sdout <- 2
normadj <- (0.5+0.5*bicor(exprs, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(exprs)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(meta$Group), pch=19, main="Outlier detection", ylab="Network Connectivity (z score)")
abline(h=-2, lty=2)
exprs = exprs[,!outliers]
meta = meta[!outliers,]

#combat
mod = model.matrix(~Group, data=meta)
batch = factor(meta$Batch)
exprs = ComBat(exprs, batch=batch, mod=mod)

#collapse rows
realGenes = !is.na(annot$ensembl_gene_id)  #
exprs = exprs[realGenes,]
annot = annot[realGenes,]

CR = collapseRows(exprs, rowGroup = annot$ensembl_gene_id, rowID = annot$affy_hugene_1_0_st_v1) 
exprs = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], annot$affy_hugene_1_0_st_v1)
annot = annot[idx,]
rownames(annot) = annot$ensembl_gene_id

#covariate regression
X = model.matrix(~Group +Sex + Age + pH + PMI + RNAdeg, meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,])))  # PMI + Sex +  Age + RNAdeg
exprs = exprs - t(to_regress)

#randomize ctl
idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA

#save
save(file="./working_data/01_unbiased/chen_normalized.RData", exprs, meta, annot)

