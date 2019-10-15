#ALTERNATE NORMALIZATION WITHOUT REGRESSION FOR SEX COVARIATE####
#IWAMOTO####

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

#02 normalization and regression, without sex####
load("./working_data/01_unbiased/iwamoto_compiled.RData")

#remove singular batches, sub quality
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
plot(z.ku, col = as.numeric(meta$Group), pch=19, main="Outlier detection", ylab="Network Connectivity (z score)")
abline(h=-2, lty=2)
exprs = exprs[,!outliers]
meta = meta[!outliers,]

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
dim(exprs)

#covariate regression
X = model.matrix(~Group+Age+pH+PMI+RNAdeg,data=meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:7]) %*% (as.matrix(beta[4:7,]))) 
exprs = exprs - t(to_regress)

idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/01_unbiased/iwamoto_normalized_sex.RData", exprs, meta, annot)
