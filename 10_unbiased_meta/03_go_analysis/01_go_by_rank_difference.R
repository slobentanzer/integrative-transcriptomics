#rank differences, GO analysis####
rm(list=ls())
options(stringsAsFactors=F)

library(nlme)

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(topGO)
library(RNeo4j)
library(biomaRt)
graph <- startGraph("http://localhost:7474/db/data/")
rank_histo <- list()
load("working_data/unbiased/meta_for_GO.RData")

get.entrez.by.score.decreasing <- function(x) {
  x <- x[order(x, decreasing = T)]
  return(as.numeric(names(x[x>cutoff])))
}

geneID2GO <- readRDS(file = "working_data/geneID2GO_all.rds")

#BETWEEN BD AND SCZ, NO SEX DIFFERENTIATION####
scz_bd_comb$rank_scz <- rank(scz_bd_comb$beta_scz, na.last = F)
scz_bd_comb$rank_bd <- rank(scz_bd_comb$beta_bd, na.last = F)
scz_bd_comb$rank_dif <- abs(scz_bd_comb$rank_scz - scz_bd_comb$rank_bd)
scz_bd_comb <- scz_bd_comb[order(scz_bd_comb$rank_dif, decreasing = T),]
hist(scz_bd_comb$rank_dif)

genes <- scz_bd_comb$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))


#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
cutoff <- quantile(allGenes, seq(0,1,.01))[96]
topGenes <- function(allScore) {
  return(allScore >= cutoff)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test 
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test 
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total 
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.bd <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.bd #inflammatory processes, bacterial defense, organophosphate

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.bd <- lapply(scoresInTerm(GOdata, res.scz.bd$GO.ID, use.names = T), get.entrez.by.score.decreasing)

#organophos
cutoff <- quantile(allGenes, seq(0,1,.01))[94]
topGenes <- function(allScore) {
  return(allScore >= cutoff)
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test 
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test 
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total 
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.bd94 <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.bd94 #inflammatory processes, bacterial defense, organophosphate


genes.scz.bd94 <- lapply(scoresInTerm(GOdata, res.scz.bd94$GO.ID, use.names = T), get.entrez.by.score.decreasing)



#BETWEEN BD AND SCZ, SCZ MINUS BD -> BD-ELEVATED####
scz_bd_comb$rank_scz <- rank(scz_bd_comb$beta_scz, na.last = F)
scz_bd_comb$rank_bd <- rank(scz_bd_comb$beta_bd, na.last = F)
scz_bd_comb$rank_dif <- scz_bd_comb$rank_scz - scz_bd_comb$rank_bd
scz_bd_comb <- scz_bd_comb[order(scz_bd_comb$rank_dif, decreasing = T),]
hist(scz_bd_comb$rank_dif)
#histogram
rank_histo$scz_bd <- scz_bd_comb$rank_dif

genes <- scz_bd_comb$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test 
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total 
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.minus.bd <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.minus.bd #neutrophils, endocytosis, phagocytosis, inflammation, complement


cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.minus.bd <- lapply(scoresInTerm(GOdata, res.scz.minus.bd$GO.ID, use.names = T), get.entrez.by.score.decreasing)

#no organophos

#BETWEEN BD AND SCZ, BD MINUS SCZ -> SCZ-ELEVATED####
scz_bd_comb$rank_scz <- rank(scz_bd_comb$beta_scz, na.last = F)
scz_bd_comb$rank_bd <- rank(scz_bd_comb$beta_bd, na.last = F)
scz_bd_comb$rank_dif <- scz_bd_comb$rank_bd - scz_bd_comb$rank_scz
scz_bd_comb <- scz_bd_comb[order(scz_bd_comb$rank_dif, decreasing = T),]
hist(scz_bd_comb$rank_dif)

genes <- scz_bd_comb$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.minus.scz <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.minus.scz #hormone, development

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.minus.scz <- lapply(scoresInTerm(GOdata, res.bd.minus.scz$GO.ID, use.names = T), get.entrez.by.score.decreasing)


#organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[99])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.minus.scz99 <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.minus.scz99 #hormone, development, organophosphate

cutoff <- quantile(allGenes, seq(0,1,.01))[99]
genes.bd.minus.scz99 <- lapply(scoresInTerm(GOdata, res.bd.minus.scz99$GO.ID, use.names = T), get.entrez.by.score.decreasing)








#HIGHER IN BD####
#BETWEEN BD AND SCZ - male - BD ELEVATED####
scz_bd_comb_m <- cbind(scz_meta_dge_m, bd_meta_dge_m)
colnames(scz_bd_comb_m) <- paste0(colnames(scz_bd_comb_m), rep(c("_scz", "_bd"), each = 5))
scz_bd_comb_m$rank_scz <- rank(scz_bd_comb_m$beta_scz, na.last = F)
scz_bd_comb_m$rank_bd <- rank(scz_bd_comb_m$beta_bd, na.last = F)
scz_bd_comb_m$rank_dif <- scz_bd_comb_m$rank_scz - scz_bd_comb_m$rank_bd
scz_bd_comb_m <- scz_bd_comb_m[order(scz_bd_comb_m$rank_dif, decreasing = T),]
hist(scz_bd_comb_m$rank_dif)
#histogram
rank_histo$scz_bd_male <- scz_bd_comb_m$rank_dif


genes <- scz_bd_comb_m$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb_m)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ, male"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.minus.bd_m <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.minus.bd_m #neutrophil, inflammation, tnf, leukocytes, phagocytosis, cytokines, IG, LPS, IL-6


#no organophos



cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.minus.bd_m <- lapply(scoresInTerm(GOdata, res.scz.minus.bd_m$GO.ID, use.names = T), get.entrez.by.score.decreasing)









#BETWEEN BD AND SCZ - female - BD ELEVATED####
scz_bd_comb_f <- cbind(scz_meta_dge_f, bd_meta_dge_f)
colnames(scz_bd_comb_f) <- paste0(colnames(scz_bd_comb_f), rep(c("_scz", "_bd"), each = 5))
scz_bd_comb_f$rank_scz <- rank(scz_bd_comb_f$beta_scz, na.last = F)
scz_bd_comb_f$rank_bd <- rank(scz_bd_comb_f$beta_bd, na.last = F)
scz_bd_comb_f$rank_dif <- scz_bd_comb_f$rank_scz - scz_bd_comb_f$rank_bd
scz_bd_comb_f <- scz_bd_comb_f[order(scz_bd_comb_f$rank_dif, decreasing = T),]
hist(scz_bd_comb_f$rank_dif)
#histogram
rank_histo$scz_bd_female <- scz_bd_comb_f$rank_dif

genes <- scz_bd_comb_f$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb_f)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ, female"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.minus.bd_f <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.minus.bd_f #phagocytosis, leukocytes, macrophages, inflammation, jak-stat, 


cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.minus.bd_f <- lapply(scoresInTerm(GOdata, res.scz.minus.bd_f$GO.ID, use.names = T), get.entrez.by.score.decreasing)




#organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[91])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "classic", ranksOf = "classic", topNodes = 100)
res.scz.minus.bd_f91 <- allRes[as.numeric(allRes$classic)<.05,]
res.scz.minus.bd_f91 #phagocytosis, leukocytes, macrophages, inflammation, jak-stat, 

#organophos only in classic, low

cutoff <- quantile(allGenes, seq(0,1,.01))[91]
genes.scz.minus.bd_f91 <- lapply(scoresInTerm(GOdata, res.scz.minus.bd_f91$GO.ID, use.names = T), get.entrez.by.score.decreasing)




#HIGHER IN SCZ####
#BETWEEN BD AND SCZ - male - SCZ ELEVATED####
scz_bd_comb_m <- cbind(scz_meta_dge_m, bd_meta_dge_m)
colnames(scz_bd_comb_m) <- paste0(colnames(scz_bd_comb_m), rep(c("_scz", "_bd"), each = 5))
scz_bd_comb_m$rank_scz <- rank(scz_bd_comb_m$beta_scz, na.last = F)
scz_bd_comb_m$rank_bd <- rank(scz_bd_comb_m$beta_bd, na.last = F)
scz_bd_comb_m$rank_dif <- scz_bd_comb_m$rank_bd - scz_bd_comb_m$rank_scz
scz_bd_comb_m <- scz_bd_comb_m[order(scz_bd_comb_m$rank_dif, decreasing = T),]
hist(scz_bd_comb_m$rank_dif)

genes <- scz_bd_comb_m$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb_m)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ, male"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.minus.scz_m <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.minus.scz_m #neurotransmitters, hormones, inflammation


cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.minus.scz_m <- lapply(scoresInTerm(GOdata, res.bd.minus.scz_m$GO.ID, use.names = T), get.entrez.by.score.decreasing)



#organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[100])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.minus.scz_m100 <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.minus.scz_m100 #neurotransmitters, hormones, inflammation



cutoff <- quantile(allGenes, seq(0,1,.01))[100]
genes.bd.minus.scz_m100 <- lapply(scoresInTerm(GOdata, res.bd.minus.scz_m100$GO.ID, use.names = T), get.entrez.by.score.decreasing)


#BETWEEN BD AND SCZ - female - SCZ ELEVATED####
scz_bd_comb_f <- cbind(scz_meta_dge_f, bd_meta_dge_f)
colnames(scz_bd_comb_f) <- paste0(colnames(scz_bd_comb_f), rep(c("_scz", "_bd"), each = 5))
scz_bd_comb_f$rank_scz <- rank(scz_bd_comb_f$beta_scz, na.last = F)
scz_bd_comb_f$rank_bd <- rank(scz_bd_comb_f$beta_bd, na.last = F)
scz_bd_comb_f$rank_dif <- scz_bd_comb_f$rank_bd - scz_bd_comb_f$rank_scz
scz_bd_comb_f <- scz_bd_comb_f[order(scz_bd_comb_f$rank_dif, decreasing = T),]
hist(scz_bd_comb_f$rank_dif)

genes <- scz_bd_comb_f$rank_dif[1:2000]
ensgs <- rownames(scz_bd_comb_f)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference BD vs SCZ, female"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.minus.scz_f <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.minus.scz_f #neuron differentiation, development, growth factor, erk, wnt, neuron fate, circadian rhythm, learning/memory

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.minus.scz_f <- lapply(scoresInTerm(GOdata, res.bd.minus.scz_f$GO.ID, use.names = T), get.entrez.by.score.decreasing)


##BETWEEN MALES AND FEMALES####
##BD####
bd_meta_comb$rank_f <- rank(bd_meta_comb$beta_f, na.last = F)
bd_meta_comb$rank_m <- rank(bd_meta_comb$beta_m, na.last = F)
bd_meta_comb$rank_dif <- abs(bd_meta_comb$rank_f - bd_meta_comb$rank_m)
bd_meta_comb <- bd_meta_comb[order(bd_meta_comb$rank_dif, decreasing = T),]
hist(bd_meta_comb$rank_dif)

genes <- bd_meta_comb$rank_dif[1:2000]
ensgs <- rownames(bd_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[1:10,]
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference between M and F in BD"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.sex <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.sex #neutrophils, fatty acids, alcohol

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.sex <- lapply(scoresInTerm(GOdata, res.bd.sex$GO.ID, use.names = T), get.entrez.by.score.decreasing)

#no organophos

##BD, FEMALE MINUS MALE -> MALE ELEVATED####
bd_meta_comb$rank_f <- rank(bd_meta_comb$beta_f, na.last = F)
bd_meta_comb$rank_m <- rank(bd_meta_comb$beta_m, na.last = F)
bd_meta_comb$rank_dif <- bd_meta_comb$rank_f - bd_meta_comb$rank_m
bd_meta_comb <- bd_meta_comb[order(bd_meta_comb$rank_dif, decreasing = T),]
hist(bd_meta_comb$rank_dif)
#histogram
rank_histo$female_male_bd <- bd_meta_comb$rank_dif

genes <- bd_meta_comb$rank_dif[1:2000]
ensgs <- rownames(bd_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[1:10,]
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference between M and F in BD"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.f.minus.m <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.f.minus.m #progesterone, alcohol, circadian rhythm, erk, inflammation, immune processes


cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.f.minus.m <- lapply(scoresInTerm(GOdata, res.bd.f.minus.m$GO.ID, use.names = T), get.entrez.by.score.decreasing)



# organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[97])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.f.minus.m97 <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.f.minus.m97 #progesterone, alcohol, circadian rhythm, erk, inflammation, immune processes

cutoff <- quantile(allGenes, seq(0,1,.01))[97]
genes.bd.f.minus.m97 <- lapply(scoresInTerm(GOdata, res.bd.f.minus.m97$GO.ID, use.names = T), get.entrez.by.score.decreasing)



##BD, MALE MINUS FEMALE -> FEMALE ELEVATED####
bd_meta_comb$rank_f <- rank(bd_meta_comb$beta_f, na.last = F)
bd_meta_comb$rank_m <- rank(bd_meta_comb$beta_m, na.last = F)
bd_meta_comb$rank_dif <- bd_meta_comb$rank_m - bd_meta_comb$rank_f
bd_meta_comb <- bd_meta_comb[order(bd_meta_comb$rank_dif, decreasing = T),]
hist(bd_meta_comb$rank_dif)

genes <- bd_meta_comb$rank_dif[1:2000]
ensgs <- rownames(bd_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[1:10,]
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference between M and F in BD"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.bd.m.minus.f <- allRes[as.numeric(allRes$weight)<.05,]
res.bd.m.minus.f #spermatogenesis, development, hormones, transmitters

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.bd.m.minus.f <- lapply(scoresInTerm(GOdata, res.bd.m.minus.f$GO.ID, use.names = T), get.entrez.by.score.decreasing)

#no organophos


#SCZ####
scz_meta_comb$rank_f <- rank(scz_meta_comb$beta_f, na.last = F)
scz_meta_comb$rank_m <- rank(scz_meta_comb$beta_m, na.last = F)
scz_meta_comb$rank_dif <- abs(scz_meta_comb$rank_f - scz_meta_comb$rank_m)
scz_meta_comb <- scz_meta_comb[order(scz_meta_comb$rank_dif, decreasing = T),]
hist(scz_meta_comb$rank_dif)

genes <- scz_meta_comb$rank_dif[1:2000]
ensgs <- rownames(scz_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference M vs F in SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.sex <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.sex #antigen processing, t cell activation, IG, cytokine, leukocytes, erk, jak-stat, neuronal apoptosis

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.sex <- lapply(scoresInTerm(GOdata, res.scz.sex$GO.ID, use.names = T), get.entrez.by.score.decreasing)

#organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[93])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.sex93 <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.sex93 #antigen processing, t cell activation, IG, cytokine, leukocytes, erk, jak-stat, neuronal apoptosis

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.sex93 <- lapply(scoresInTerm(GOdata, res.scz.sex93$GO.ID, use.names = T), get.entrez.by.score.decreasing)




#SCZ, FEMALE MINUS MALE -> MALE ELEVATED####
scz_meta_comb$rank_f <- rank(scz_meta_comb$beta_f, na.last = F)
scz_meta_comb$rank_m <- rank(scz_meta_comb$beta_m, na.last = F)
scz_meta_comb$rank_dif <- scz_meta_comb$rank_f - scz_meta_comb$rank_m
scz_meta_comb <- scz_meta_comb[order(scz_meta_comb$rank_dif, decreasing = T),]
hist(scz_meta_comb$rank_dif)
#histogram
rank_histo$female_male_scz <- scz_meta_comb$rank_dif

genes <- scz_meta_comb$rank_dif[1:2000]
ensgs <- rownames(scz_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference M vs F in SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.f.minus.m <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.f.minus.m #endocytosis, viral entry, t cells, cytokines, antigen, erk, interferon, jak-stat

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.f.minus.m <- lapply(scoresInTerm(GOdata, res.scz.f.minus.m$GO.ID, use.names = T), get.entrez.by.score.decreasing)



#organophos
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[97])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "classic", ranksOf = "classic", topNodes = 100)
res.scz.f.minus.m97 <- allRes[as.numeric(allRes$classic)<.05,]
res.scz.f.minus.m97 #endocytosis, viral entry, t cells, cytokines, antigen, erk, interferon, jak-stat

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.f.minus.m97 <- lapply(scoresInTerm(GOdata, res.scz.f.minus.m97$GO.ID, use.names = T), get.entrez.by.score.decreasing)



#SCZ, MALE MINUS FEMALE -> FEMALE ELEVATED####
scz_meta_comb$rank_f <- rank(scz_meta_comb$beta_f, na.last = F)
scz_meta_comb$rank_m <- rank(scz_meta_comb$beta_m, na.last = F)
scz_meta_comb$rank_dif <- scz_meta_comb$rank_m - scz_meta_comb$rank_f
scz_meta_comb <- scz_meta_comb[order(scz_meta_comb$rank_dif, decreasing = T),]
hist(scz_meta_comb$rank_dif)

genes <- scz_meta_comb$rank_dif[1:2000]
ensgs <- rownames(scz_meta_comb)[1:2000]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensg2entrez <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                     filters = "ensembl_gene_id",
                     values = ensgs,
                     mart = mart)
names(genes) <- ensg2entrez$entrezgene[match(ensgs, ensg2entrez$ensembl_gene_id)]
which(is.na(names(genes)))

geneID2GO <- readRDS(file = "../250118_miRNA_seq/data/geneID2GO.rds")

#ontology
ontology <- "BP"
#description
description <- "Difference M vs F in SCZ"
#allGenes list
allGenes <- genes
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
topGenes <- function(allScore) {
  return(allScore >= quantile(allGenes, seq(0,1,.01))[96])
}
#prune GO terms
nodeSize <- 10

GOdata <- new("topGOdata",
              description = description,
              ontology = ontology,
              allGenes = allGenes,
              geneSel = topGenes,
              annot = annFUN.gene2GO,
              nodeSize = nodeSize,
              gene2GO = geneID2GO)

## >fisher test ###
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ###
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

## >total ###
allRes <- GenTable(GOdata, classic = resultFisher, 
                   #KS = resultKS, 
                   weight = resultWeight, #elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
res.scz.m.minus.f <- allRes[as.numeric(allRes$weight)<.05,]
res.scz.m.minus.f #almost nothing

cutoff <- quantile(allGenes, seq(0,1,.01))[96]
genes.scz.m.minus.f <- lapply(scoresInTerm(GOdata, res.scz.m.minus.f$GO.ID, use.names = T), get.entrez.by.score.decreasing)


#no organophos


#save####
save(res.scz.bd, #inflammatory processes, bacterial defense, development, cytokines, organophosphate
     res.scz.bd94,
     res.scz.minus.bd,
     res.bd.minus.scz, #differentiation and development, hormones, wnt
     res.bd.minus.scz99,
     res.scz.minus.bd_m,
     res.scz.minus.bd_f, #phagocytosis, leukocytes, macrophages, inflammation, pregnancy, LPS, jak-stat, bacterial response
     res.scz.minus.bd_f91,
     res.bd.minus.scz_m, #neurotransmitters, hormones, inflammation, telomere maintenance, fertilization
     res.bd.minus.scz_m100,
     res.bd.minus.scz_f,
     res.bd.sex,
     res.bd.f.minus.m,
     res.bd.f.minus.m97,
     res.bd.m.minus.f,
     res.scz.sex,
     res.scz.sex93,
     res.scz.f.minus.m, #endocytosis, viral entry, t cells, cytokines, antigen, erk, enuron death, inflammation, LPS, aging, oxidative stress
     res.scz.f.minus.m97,
     res.scz.m.minus.f, file = "working_data/01_unbiased/GO_tables.RData")

save(genes.scz.bd, #inflammatory processes, bacterial defense, development, cytokines, organophosphate
     genes.scz.bd94,
     genes.scz.minus.bd,
     genes.bd.minus.scz, #differentiation and development, hormones, wnt
     genes.bd.minus.scz99,
     genes.scz.minus.bd_m,
     genes.scz.minus.bd_f, #phagocytosis, leukocytes, macrophages, inflammation, pregnancy, LPS, jak-stat, bacterial genesponse
     genes.scz.minus.bd_f91,
     genes.bd.minus.scz_m, #neurotransmitters, hormones, inflammation, telomere maintenance, fertilization
     genes.bd.minus.scz_m100,
     genes.bd.minus.scz_f,
     genes.bd.sex,
     genes.bd.f.minus.m,
     genes.bd.f.minus.m97,
     genes.bd.m.minus.f,
     genes.scz.sex,
     genes.scz.sex93,
     genes.scz.f.minus.m, #endocytosis, viral entry, t cells, cytokines, antigen, erk, enuron death, inflammation, LPS, aging, oxidative stgeness
     genes.scz.f.minus.m97,
     genes.scz.m.minus.f, file = "working_data/01_unbiased/GO_genes.RData")


# RESULTS####
load("working_data/01_unbiased/GO_genes.RData")
load("working_data/01_unbiased/GO_tables.RData")

#DISEASE WITHOUT SEX DIFFERENTIATION
#general scz vs bd, organophos
res.scz.bd #inflammatory processes, bacterial defense, development, cytokines, organophosphate
res.scz.bd94

#elevated in BD, no o
res.scz.minus.bd #neutrophils, endocytosis, phagocytosis, inflammation, complement
#elevated in SCZ, organophos top 2%
res.bd.minus.scz #differentiation and development, hormones, wnt
res.bd.minus.scz99 #glucose, metabolism (starvation), wnt, organophosphate, differentiation, fertilization

#intersection
res.scz.minus.bd[res.scz.minus.bd$GO.ID %in% res.bd.minus.scz$GO.ID,] #nothing

#DISEASE WITH SEX DIFFERENTIATION
#per sex elevated in BD
#male, no o
res.scz.minus.bd_m #neutrophil, inflammation, tnf, leukocytes, phagocytosis, cytokines, IG, LPS, IL-6, NFKB
#female, organophos only in classic, low
res.scz.minus.bd_f #phagocytosis, leukocytes, macrophages, inflammation, pregnancy, LPS, jak-stat, bacterial response
res.scz.minus.bd_f91 #lots of inflammation, organophosphates

#intersection
res.scz.minus.bd_m[res.scz.minus.bd_m$GO.ID %in% res.scz.minus.bd_f$GO.ID,] #inflammation, tnf, leukocytes, bacterial response

#not in non-sex-differentiated
res.scz.minus.bd_m[!res.scz.minus.bd_m$GO.ID %in% res.scz.minus.bd$GO.ID,] #
res.scz.minus.bd_f[!res.scz.minus.bd_f$GO.ID %in% res.scz.minus.bd$GO.ID,] #

#per sex elevated in SCZ
#male, organophos top 1%
res.bd.minus.scz_m #neurotransmitters, hormones, inflammation, telomere maintenance, fertilization
res.bd.minus.scz_m100 #hormones, differentiation, telomere maintenance, organophosphates, viral response
#female, no organophos
res.bd.minus.scz_f #neuron differentiation, development, growth factor, erk cascade, wnt, neuron fate, starvation, 
# circadian rhythm, hypoxia, learning/memory

#intersection
res.bd.minus.scz_m[res.bd.minus.scz_m$GO.ID %in% res.bd.minus.scz_f$GO.ID,] #nothing

#not in non-sex-differentiated
res.bd.minus.scz_m[!res.bd.minus.scz_m$GO.ID %in% res.bd.minus.scz$GO.ID,] #neurotransmitters, hormones, inflammation, nutrients, 
# telomere maintenance
res.bd.minus.scz_f[!res.bd.minus.scz_f$GO.ID %in% res.bd.minus.scz$GO.ID,] #neuron differentiation, development, growth factor, 
# erk cascade, neuron fate, circadian rhythm, hypoxia, learning/memory

#DIFFERENCE BETWEEN SEXES
#general sex diff bd, no o
res.bd.sex #neutrophils, fatty acids, alcohol, cognition, tnf
#elevated in males, organophos
res.bd.f.minus.m #progesterone, alcohol, circadian rhythm, cell fate, oligo differentiation, cell maturation, 
# erk cascade, inflammation, gliogenesis, immune processes
res.scz.f.minus.m #endocytosis, viral entry, t cells, cytokines, antigen, erk, neuron death, inflammation, LPS, aging, oxidative stress
res.bd.f.minus.m97 #alcohol, circadian rhythm, cell fate, erk cascade, fatty acids, macrophages, organophosphates, 
# glial differentiation, immune system, inflammation, oligo differentiation
#elevated in females, no o
res.bd.m.minus.f #spermatogenesis, development, fertilization, female gamete generation, hormones, transmitters,
# developmental maturation, neuron development

#intersection
res.bd.f.minus.m[res.bd.f.minus.m$GO.ID %in% res.bd.m.minus.f$GO.ID,] #almost nothing
res.bd.f.minus.m[res.bd.f.minus.m$GO.ID %in% res.scz.f.minus.m$GO.ID,] #regulated exocytosis
res.bd.m.minus.f[res.bd.m.minus.f$GO.ID %in% res.scz.m.minus.f$GO.ID,] #nothing

#general sex diff scz, organophos
res.scz.sex #antigen processing, t cell activation, IG, cytokine, leukocytes, erk cascade, jak-stat, inf, 
# inflammatory response, neuronal apoptosis
res.scz.sex93 #IG, jak-stat, organophosphates, immune response, proliferation, inflammation, cytokines, glial cell development
#elevated in males, organophosphates only in classic, low
res.scz.f.minus.m #endocytosis, viral entry, t cells, cytokines, antigen, erk, neuron death, inflammation, LPS, aging, oxidative stress
res.scz.f.minus.m97 #viral entry, cytokine, immune response, organophosphate
#elevated in females, no o
res.scz.m.minus.f #almost nothing

#intersection
res.scz.f.minus.m[res.scz.f.minus.m$GO.ID %in% res.scz.m.minus.f$GO.ID,] #nothing

#get individual genes per GO term####
result_list <- list(res.scz.bd, #inflammatory processes, bacterial defense, development, cytokines, organophosphate
                    res.scz.bd94,
                    res.scz.minus.bd,
                    res.bd.minus.scz, #differentiation and development, hormones, wnt
                    res.bd.minus.scz99,
                    res.scz.minus.bd_m,
                    res.scz.minus.bd_f, #phagocytosis, leukocytes, macrophages, inflammation, pregnancy, LPS, jak-stat, bacterial response
                    res.scz.minus.bd_f91,
                    res.bd.minus.scz_m, #neurotransmitters, hormones, inflammation, telomere maintenance, fertilization
                    res.bd.minus.scz_m100,
                    res.bd.minus.scz_f,
                    res.bd.sex,
                    res.bd.f.minus.m,
                    res.bd.f.minus.m97,
                    res.bd.m.minus.f,
                    res.scz.sex,
                    res.scz.sex93,
                    res.scz.f.minus.m, #endocytosis, viral entry, t cells, cytokines, antigen, erk, enuron death, inflammation, LPS, aging, oxidative stress
                    res.scz.f.minus.m97,
                    res.scz.m.minus.f)

gene_list <- list(genes.scz.bd, #inflammatory processes, bacterial defense, development, cytokines, organophosphate
                  genes.scz.bd94,
                  genes.scz.minus.bd,
                  genes.bd.minus.scz, #differentiation and development, hormones, wnt
                  genes.bd.minus.scz99,
                  genes.scz.minus.bd_m,
                  genes.scz.minus.bd_f, #phagocytosis, leukocytes, macrophages, inflammation, pregnancy, LPS, jak-stat, bacterial genesponse
                  genes.scz.minus.bd_f91,
                  genes.bd.minus.scz_m, #neurotransmitters, hormones, inflammation, telomere maintenance, fertilization
                  genes.bd.minus.scz_m100,
                  genes.bd.minus.scz_f,
                  genes.bd.sex,
                  genes.bd.f.minus.m,
                  genes.bd.f.minus.m97,
                  genes.bd.m.minus.f,
                  genes.scz.sex,
                  genes.scz.sex93,
                  genes.scz.f.minus.m, #endocytosis, viral entry, t cells, cytokines, antigen, erk, enuron death, inflammation, LPS, aging, oxidative stgeness
                  genes.scz.f.minus.m97,
                  genes.scz.m.minus.f)

condition_names <- c("SCZ-BD", "SCZ-BD-94", "BD", "SCZ", "SCZ-99", "BD-m", "BD-f", "BD-f-91", "SCZ-m", "SCZ-m-100", "SCZ-f", 
                     "SEX-BD", "M-BD", "M-BD-97", "F-BD", "SEX-SCZ", "SEX-SCZ-93", "M-SCZ", "M-SCZ-97", "F-SCZ")

genes_and_terms_per_condition <- vector(mode = "list", length = length(result_list))
for(j in 1:length(result_list)) {
  res <- result_list[[j]]
  genes <- gene_list[[j]]
  n <- length(res$GO.ID)
  genes.per.term <- vector(mode = "list", length = n)
  for (i in 1:n) {
    id <- res$GO.ID[i]
    term <- res$Term[i]
    entrez <- genes[id][[1]]
    if (length(entrez) > 1)
      gene_frame <-
      cypher(graph,
             "MATCH (g:GENE) WHERE g.id IN {ids} RETURN g.name, g.id",
             ids = entrez)
    else
      gene_frame <-
      cypher(graph,
             paste0("MATCH (g:GENE) WHERE g.id = ", entrez, " RETURN g.name, g.id"))
    gene_frame <- gene_frame[match(entrez, gene_frame$g.id), ]
    rownames(gene_frame) <- seq(nrow(gene_frame))
    genes.per.term[[i]] <- gene_frame
    names(genes.per.term)[i] <- paste(id, term)
  }
  genes_and_terms_per_condition[[j]] <- genes.per.term
  names(genes_and_terms_per_condition)[j] <- condition_names[j]
}

idx <- grep("BD", names(genes_and_terms_per_condition))
for(i in seq_along(genes_and_terms_per_condition[idx])) {
  print(names(genes_and_terms_per_condition[idx])[i])
  list <- genes_and_terms_per_condition[idx][[i]]
  print(list)
  readline(prompt = "Press [ENTER]")
}

#find genes of interest####
load("../00_data/ecs_tables.RData")
mir_targets <- readRDS("../00_data/mir_sequencing_top_gene_targets_weighted.rds")
goi_entrez <- as.numeric(c(cholinergic_core_genes$entrez[cholinergic_core_genes$group != "ffar"], mir_targets$g.entrez[mir_targets$top]))

idx <- grep("BD", names(genes_and_terms_per_condition))
for(i in seq_along(genes_and_terms_per_condition[idx])) {
  print(names(genes_and_terms_per_condition[idx])[i])
  list <- genes_and_terms_per_condition[idx][[i]]
  
  print(list[which(unlist(lapply(list, function(x) any(cholinergic_core_genes$entrez[
    cholinergic_core_genes$group %in% c("nt", "nt_receptor")] %in% x[,2]))))])
  readline(prompt = "Press [ENTER]")
}

#write output
write.table(res.scz.bd, "out/sheet1.csv", row.names = F)
write.table(res.scz.minus.bd, "out/sheet2.csv", row.names = F)
write.table(res.bd.minus.scz, "out/sheet3.csv", row.names = F)
write.table(res.scz.minus.bd_m, "out/sheet4.csv", row.names = F)
write.table(res.scz.minus.bd_f, "out/sheet5.csv", row.names = F)
write.table(res.bd.minus.scz_m, "out/sheet6.csv", row.names = F)
write.table(res.bd.minus.scz_f, "out/sheet7.csv", row.names = F)
write.table(res.bd.sex, "out/sheet8.csv", row.names = F)
write.table(res.bd.f.minus.m, "out/sheet9.csv", row.names = F)
write.table(res.bd.m.minus.f, "out/sheet10.csv", row.names = F)
write.table(res.scz.sex, "out/sheet11.csv", row.names = F)
write.table(res.scz.f.minus.m, "out/sheet12.csv", row.names = F)
write.table(res.scz.m.minus.f, "out/sheet13.csv", row.names = F)

#histo plot####
saveRDS(rank_histo, file = "working_data/GO_rank_histograms.rds")
rank_histo <- readRDS(file = "working_data/GO_rank_histograms.rds")
rank_histo_frame <- data.frame()
for(i in 1:length(rank_histo)) {
  temp <- rank_histo[[i]]
  df <- data.frame(set = names(rank_histo)[i], values = temp)
  rank_histo_frame <- rbind(rank_histo_frame, df)
}

rank_histo_frame$set <- factor(rank_histo_frame$set, levels = c("scz_bd", "scz_bd_female", "scz_bd_male", "female_male_bd", "female_male_scz"))
levels(rank_histo_frame$set) <- c("SCZ vs BD", "SCZ vs BD, females", "SCZ vs BD, males", "F vs M, BD", "F vs M, SCZ")

ggplot(rank_histo_frame, aes(values, color = set)) + geom_freqpoly(binwidth = 100, alpha = .6) +
  xlab("Spearman's rank difference") + ylab("Frequency") + scale_color_brewer(palette = "Dark2", name = "Comparison")
ggsave("img/dge_rank_histogram.pdf", height = 4, width = 10)
