#GENE TARGET PREDICTION OF CNTF-INDUCED MIRS####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(ggplot2)
library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")


mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")

long_term <- rownames(mir.matrix)[rowSums(mir.matrix[, 1:8]) != 0]

# long_term <- unlist(lapply(strsplit(long_term, "-"), function(x) paste(x[1:3], collapse = "-")))
# long_term <- gsub("r", "R", long_term)
# long_term <- unique(long_term)


targets <- cypher(graph, "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE) WHERE m.name in {names} 
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
                  names = long_term)
targets <- targets[order(targets$m.name, targets$g.name),]
targets <- targets[!duplicated(paste0(targets$m.name, targets$g.ensg)),] #remove miRs with multiple ancestors
head(targets)

#GO STUFF####
#sort by occurences?
top_genes <- dplyr::count(targets, g.name, sort = T)
top_genes <- aggregate(r.rating ~ g.ensg, data = targets, FUN = "sum")
# top_genes <- dplyr::count(targets[targets$m.name %in% targets.count$m.name[grep("day", targets.count$condition)],], g.name, sort = T)# only long term
# top_genes <- dplyr::count(targets[targets$m.name %in% targets.count$m.name[grep("min", targets.count$condition)],], g.name, sort = T)# only short term
top_genes$g.name <- targets$g.name[match(top_genes$g.ensg, targets$g.ensg)]
top_genes <- top_genes[order(top_genes$r.rating, decreasing = T),]

## >>unbiased GO analysis - all results####
library(topGO)
geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")

# BPTerms <- ls(GOBPTerm)
# head(BPTerms)
#ontology
ontology <- "BP"
#description
description <- "miR targets in LA2/LA5 during CNTF differentiation (long term)"
#allGenes list
allGenes <- top_genes$r.rating[1:2000]
names(allGenes) <- top_genes$g.ensg[1:2000]
#geneSelection function
hist(allGenes)
quantile(allGenes, seq(0,1,.01))
x <- quantile(allGenes, seq(0,1,.01))[91]
topGenes <- function(allScore) {
  return(allScore >= x)
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

## >fisher test ####
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

# >weight test ####
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)
resultWeight

## >total ####
allRes <- GenTable(GOdata, classic = resultFisher, 
                   # KS = resultKS, 
                   weight = resultWeight, 
                   # elim = resultElim,
                   orderBy = "weight", ranksOf = "classic", topNodes = 100)
# allRes <- allRes[order(allRes$Annotated, decreasing = T),]
allRes[allRes$weight<.05,]

write.table(allRes[allRes$weight<.05,], file = "./out/top_mir_targets_GO_enrich.csv",
            quote = T, row.names = F)
