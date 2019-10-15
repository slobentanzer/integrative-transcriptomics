#ALTERNATE NORMALIZATION WITHOUT REGRESSION FOR SEX COVARIATE####
#LANZ####

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
load(file="./working_data/01_unbiased/lanz_compiled.RData")

#Remove singular batches
table(meta$Batch)

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

#sub quality
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
X = model.matrix(~Group+Age+PMI+pH+Race+RIN+RNAdeg, data=meta)
Y = exprs
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,5:10]) %*% (as.matrix(beta[5:10,])))  
exprs = exprs - t(to_regress)




idx=  which(meta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
meta$Group.BD = meta$Group.SCZ = meta$Group
meta$Group.SCZ[c(which(meta$Group=="BD"), idx[!subset.scz])] = NA
meta$Group.BD[c(which(meta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/01_unbiased/lanz_normalized_sex.RData", exprs, meta, annot)

  
  