#LOAD IWAMOTO DATASET, NORMALIZE, BATCH CORRECT, REGRESS COVARIATES####
#GSE12649####

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
meta = read.csv(file = "./raw_data/datasets/Iwamoto_GSE12649/Iwa_GSE12649_datMeta.csv", head =
                     T)
rownames(meta) = meta$GSM
meta$Sex = as.factor(meta$Sex)
meta$Race = as.factor(meta$Race)
meta$Group = meta$Profile
meta$Group[meta$Group == "Bipolar"] = "BD"
meta$Group[meta$Group == "Schizophrenia"] = "SCZ"
meta$Group[meta$Group == "Control"] = "CTL"
meta$Group[meta$Group == "Depression"] = "MDD"
meta$Group = factor(meta$Group, levels = c("CTL", "BD", "SCZ"))

meta$RD = as.factor(meta$RD)
meta = meta[!is.na(meta$Group), ]
meta$Study = "SCZ.BD.Iwa"

#exprs####
data.affy = ReadAffy(celfile.path = "./raw_data/datasets/Iwamoto_GSE12649/",
                     filenames = paste(meta$GSM, ".CEL.gz", sep = ""))
exprs = affy::rma(data.affy, background = T, normalize = T)
exprs = exprs(exprs)
colnames(exprs) <-
  unlist(lapply(strsplit(colnames(exprs), ".", fixed = T), "[", 1))

#meta match
idx = match(colnames(exprs), rownames(meta))
meta = meta[idx, ]

RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
meta$RNAdeg = RNAdeg$slope

sd = protocolData(data.affy)$ScanDate
sd = substring(sd, 1, nchar(sd) - 9)
meta$Batch = as.factor(sd)

#annotation (gandal)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host =
                    "feb2014.archive.ensembl.org")
identifier <- "affy_hg_u133a"
getinfo <-
  c(
    "affy_hg_u133a",
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
idx = match(rownames(exprs), geneDat$affy_hg_u133a)
annot = cbind(rownames(exprs), geneDat[idx, ])

save(file = "./working_data/01_unbiased/iwamoto_compiled.RData", exprs, meta, annot)

#02 normalization and regression####
load("./working_data/01_unbiased/iwamoto_compiled.RData")

#batch control, sub quality (gandal)
table(meta$Batch)
to_keep = !(meta$Batch == "07/09/03") & (meta$RNAdeg < 3.8) & (meta$RNAdeg > 2) & !is.na(meta$Group) & (meta$Filename != "77-HGU133A")
table(to_keep)
meta = meta[to_keep,]
exprs = exprs[,to_keep]

#remove outliers (gandal)
sdout <- 2
normadj <- (0.5+0.5*bicor(exprs, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(exprs)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(meta$Group), pch=19, main="Outlier detection", ylab="Euclidean distance (z score)")
abline(h=-2, lty=2)
exprs = exprs[,!outliers]
meta = meta[!outliers,]

meta$Batch=factor(meta$Batch)
model = model.matrix(~Group+Sex+Age+pH+PMI+RNAdeg+Batch,data=meta)

#combat
mod = model.matrix(~Group, data=meta)
batch = factor(meta$Batch)
exprs.combat = ComBat(exprs, batch=batch, mod=mod)
exprs = exprs.combat

#collapse rows
realGenes = !is.na(annot$ensembl_gene_id)  #
table(realGenes)
exprs = exprs[realGenes,]
annot = annot[realGenes,]

CR = collapseRows(exprs, rowGroup = annot$ensembl_gene_id, rowID = annot$affy_hg_u133a) 
exprs = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], annot$affy_hg_u133a)
annot = annot[idx,]
rownames(annot) = annot$ensembl_gene_id

#covariate regression
X = model.matrix(~Group+Sex+Age+pH+PMI+RNAdeg,data=meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,]))) 
exprs = exprs - t(to_regress)

#randomize ctl
idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/01_unbiased/iwamoto_normalized.RData", exprs, meta, annot)
