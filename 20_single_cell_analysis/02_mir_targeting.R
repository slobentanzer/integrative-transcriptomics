#MIR TARGETING SINGLE CELL SEQ DATA####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

genes <- cypher(graph, "MATCH (g:GENE {species: 'HSA'}) RETURN g.name")$g.name

darmanis <- readRDS("./raw_data/darmanis_chol_exprs.rds")
habib <- readRDS("./raw_data/habib_chol_exprs.rds")
tasic <- readRDS("./raw_data/tasic_chol_exprs.rds")
zeisel <- readRDS("./raw_data/zeisel_chol_exprs.rds")

#remove marker genes
marker <- rownames(darmanis)[36:42]
rowMeans(darmanis[1:35,1:9])
darmanis <- darmanis[1:35,]
rowMeans(habib[8:42,1:9])
habib <- habib[!habib$name %in% marker,]
rowMeans(tasic[!rownames(tasic) %in% marker,1:4])
tasic <- tasic[!tasic$name %in% marker,]
rowMeans(zeisel[!rownames(zeisel) %in% marker,1:4])
zeisel <- zeisel[!zeisel$name %in% marker,]

#more than one per row
darmanis <- darmanis[colSums(apply(darmanis[,1:9]>0, 1, as.numeric))>1,]
habib <- habib[colSums(apply(habib[,1:9]>0, 1, as.numeric))>1,]
tasic <- tasic[colSums(apply(tasic[,1:4]>0, 1, as.numeric))>1,]
zeisel <- zeisel[colSums(apply(zeisel[,1:4]>0, 1, as.numeric))>1,]

human <- rownames(darmanis)

mouse <- unique(c(habib$name, tasic$name, zeisel$name))

#human genes####
if(!file.exists("./working_data/single_cell/human_permutation.rds")) {
#mir targeting
targets <- cypher(graph, "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE g.name in {names} AND r.rating > 5 
                  RETURN m.name, r.rating, g.name", 
                  names = human)

human.count <- data.frame(dplyr::count(targets, m.name, sort = T))
head(human.count, 100)

#permutation
human.count$p.val <- NA
} else {
  human.count <-
    readRDS(file = "./working_data/single_cell/human_permutation.rds")
}
if(ncol(human.count)<10000) {
  for (i in (ncol(human.count) + 1):10000) {
    print(i)
    perm_genes <- sample(genes, length(human))
    targets <- cypher(
      graph,
      "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE)
      WHERE g.name in {names} AND r.rating > 5
      RETURN m.name, r.rating, g.name",
      names = perm_genes
    )
    targets.count <- data.frame(dplyr::count(targets, m.name))
    idx <- match(human.count$m.name, targets.count$m.name)
    human.count[, i] <- targets.count$n[idx]
  }
  
  human.count[is.na(human.count)] <- 0
  for (mir in rownames(human.count)) {
    null <-
      unlist(human.count[rownames(human.count) == mir, 3:ncol(human.count)])
    null <- sort(null)
    null <-
      cbind(null, data.frame(prob = 1 - abs(2 * seq(1:length(
        null
      )) / length(null) - 1)))
    count <- human.count$n[rownames(human.count) == mir]
    human.count$p.val[rownames(human.count) == mir] <-
      null$prob[findInterval(count, null$null)]
  }
  
  human.count <- human.count[order(human.count$n, decreasing = T), ]
  rownames(human.count) <- seq(nrow(human.count))
  saveRDS(human.count, file = "./working_data/single_cell/human_permutation.rds")
}
human.count[human.count$p.val<.05,c("m.name", "n", "p.val")]

#mouse genes####
if(!file.exists("./working_data/single_cell/mouse_permutation.rds")) {
  #mir targeting
  targets <- cypher(
    graph,
    "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE)
    WHERE g.name in {names} AND r.rating > 5 
    RETURN m.name, m.primate, r.rating, g.name",
    names = mouse
  )
  
  mouse.count <- data.frame(dplyr::count(targets[!targets$m.primate,], m.name, sort = T))
  head(mouse.count, 100)
  
  #permutation
  mouse.count$p.val <- NA
} else {
  mouse.count <-
    readRDS(file = "./working_data/single_cell/mouse_permutation.rds")
}
if(ncol(mouse.count)<10000) {
  for (i in (ncol(mouse.count) + 1):10000) {
    print(i)
    perm_genes <- sample(genes, length(mouse))
    targets <- cypher(
      graph,
      "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE)
      WHERE g.name in {names} AND r.rating > 5
      RETURN m.name, m.primate, r.rating, g.name",
      names = perm_genes
    )
    targets.count <-
      data.frame(dplyr::count(targets[!targets$m.primate, ], m.name, sort = T))
    idx <- match(mouse.count$m.name, targets.count$m.name)
    mouse.count[, i] <- targets.count$n[idx]
  }
  
  mouse.count[is.na(mouse.count)] <- 0
  for (mir in rownames(mouse.count)) {
    null <-
      unlist(mouse.count[rownames(mouse.count) == mir, 4:ncol(mouse.count)])
    null <- sort(null)
    null <-
      cbind(null, data.frame(prob = 1 - abs(2 * seq(1:length(
        null
      )) / length(null) - 1)))
    count <- mouse.count$n[rownames(mouse.count) == mir]
    mouse.count$p.val[rownames(mouse.count) == mir] <-
      null$prob[findInterval(count, null$null)]
  }
  
  mouse.count <- mouse.count[order(mouse.count$n, decreasing = T), ]
  rownames(mouse.count) <- seq(nrow(mouse.count))
  saveRDS(mouse.count, file = "./working_data/single_cell/mouse_permutation.rds")
}
mouse.count[mouse.count$p.val<.05,c("m.name", "n", "p.val")]

primate_specific_mirs <- human.count[human.count$p.val<.05 & 
              !human.count$m.name %in% mouse.count$m.name,
            c("m.name", "n", "p.val")]
non_ps_mirs <- human.count[human.count$p.val<.05 & 
              human.count$m.name %in% mouse.count$m.name[mouse.count$p.val<.05],
            c("m.name", "n", "p.val")]

primate_specific_mirs <- primate_specific_mirs[order(primate_specific_mirs$n, decreasing = T),]
primate_specific_mirs
nrow(primate_specific_mirs)
primate_specific_mirs[primate_specific_mirs$p.val<.001,]
nrow(primate_specific_mirs[primate_specific_mirs$p.val<.001,])
non_ps_mirs <- non_ps_mirs[order(non_ps_mirs$n, decreasing = T),]
non_ps_mirs
nrow(non_ps_mirs)
non_ps_mirs[non_ps_mirs$p.val<.001,]
nrow(non_ps_mirs[non_ps_mirs$p.val<.001,])

saveRDS(primate_specific_mirs, file = "./working_data/single_cell/primate_specific_mirs.rds")
saveRDS(non_ps_mirs, file = "./working_data/single_cell/non_primate_specific_mirs.rds")
