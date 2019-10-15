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
meta = read.csv("./raw_data/datasets/Narayan_GSE21138/Narayan_GSE21138_datMeta.csv")
rownames(meta) = meta$GE

meta$Sex = as.factor(meta$Sex)
meta$Group = factor(meta$Group, levels = c("CTL", "SCZ"))
meta$Age = as.numeric(meta$Age)
meta$PMI = as.numeric(meta$PMI)
meta$pH = as.numeric(meta$pH)
meta$Duration = as.numeric(meta$Duration)
meta$Study = "Narayan_SCZ"
meta$Subject = paste("VictorianBrainBank",
                        meta$Group,
                        meta$Age,
                        meta$Sex,
                        meta$PMI,
                        sep = "_")

#exprs####
data.affy = ReadAffy(celfile.path = "./raw_data/datasets/Narayan_GSE21138/")

RNAdeg = AffyRNAdeg(data.affy)
plotAffyRNAdeg(RNAdeg)
meta$RNAdeg = RNAdeg$slope
sd = protocolData(data.affy)$ScanDate
sd = substring(sd, 1, nchar(sd) - 9)
meta$Batch = as.factor(sd)

#rma
exprs = rma(
  data.affy,
  background = T,
  normalize = T,
  verbose = T
)
exprs = exprs(exprs)
colnames(exprs) = substr(colnames(exprs), 1, 9)
idx = match(colnames(exprs), rownames(meta))
meta = meta[idx, ]

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

save(file="./working_data/01_unbiased/narayan_compiled.RData", exprs, meta, annot)


#02 normalization and regression####
load("./working_data/01_unbiased/narayan_compiled.RData")

#batch control, sub quality
table(meta$Batch)
to_keep = (meta$Batch != "06/14/06") & (meta$RNAdeg < 5.5)
meta = meta[to_keep, ]; exprs = exprs[,to_keep]
meta$Batch = factor(meta$Batch)

#remove outliers (gandal)
sdout <- 2
normadj <- (0.5+0.5*bicor(exprs, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(exprs)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(meta$Group), pch=19, main="Outlier detection", ylab="Connectivity (Z)")
abline(h=-2, lty=2)
exprs = exprs[,!outliers]
meta = meta[!outliers,]

#combat
mod = model.matrix(~Group, data=meta)
batch = as.factor(meta$Batch)
rownames(exprs)
exprs.combat = ComBat(exprs, batch=batch, mod=mod)
rownames(exprs.combat)
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
X = model.matrix(~Group+Sex+Age+PMI+pH+RNAdeg, data=meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))
exprs = exprs - t(to_regress)

save(file = "./working_data/01_unbiased/narayan_normalized.RData", exprs, meta, annot)
