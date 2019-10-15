#DME MIR TARGETS GO ANALYSIS####
home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

## >>unbiased GO analysis - all results####
library(topGO)

get.GO <- function(scores, ensgs, quant = .9){
  geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
  
  #ontology
  ontology <- "BP"
  # #description
  # description <- "miR targets in LA-N-2 during CNTF differentiation (all conditions), miRExpress and bowtie/bwa intersect"
  #allGenes list
  allGenes <- scores
  names(allGenes) <- ensgs
  #geneSelection function
  hist(allGenes)
  quantile(allGenes, seq(0, 1, .01))
  x <- quantile(allGenes, seq(0, 1, quant))
  topGenes <- function(allScore) {
    return(allScore >= x)
  }
  #prune GO terms
  nodeSize <- 10
  
  GOdata <- new(
    "topGOdata",
    # description = description,
    ontology = ontology,
    allGenes = allGenes,
    geneSel = topGenes,
    annot = annFUN.gene2GO,
    nodeSize = nodeSize,
    gene2GO = geneID2GO
  )
  
  ## >fisher test ####
  test.stat <-
    new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  resultFisher
  
  # >weight test ####
  test.stat <-
    new(
      "weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test",
      sigRatio = "ratio"
    )
  resultWeight <- getSigGroups(GOdata, test.stat)
  resultWeight
  
  ## >total ####
  allRes <- GenTable(
    GOdata,
    classic = resultFisher,
    weight = resultWeight,
    orderBy = "weight",
    ranksOf = "classic",
    topNodes = 100
  )
  res <- allRes[allRes$weight < .05, ]
  
  return(list(res, GOdata))
}