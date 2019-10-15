#ALTERNATE NORMALIZATION WITHOUT REGRESSION FOR SEX COVARIATE####
#CHEN####

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

#sub quality
to_keep = meta$Group!="MDD" & meta$PMI < 70
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
X = model.matrix(~Group + Age + pH + PMI + RNAdeg, meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:7]) %*% (as.matrix(beta[4:7,])))  # PMI +  Age + RNAdeg
exprs = exprs - t(to_regress)



idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA


save(file="./working_data/01_unbiased/chen_normalized_sex.RData", exprs, meta, annot)

