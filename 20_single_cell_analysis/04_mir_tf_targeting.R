#MIR TARGETING TFs PREDICTED FROM SINGLE CELL SEQ DATA####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

primate_specific_mirs <- readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <- readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
sig_tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")

genes <- cypher(graph, "MATCH (g:GENE {species: 'HSA', tf: TRUE}) RETURN g.name")$g.name

#human genes####
if(!file.exists("./working_data/single_cell/mir_tf_permutation.rds")) {
#mir targeting
targets <- cypher(graph, "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name in {mirs} AND g.name in {genes} AND r.rating > 5 
                  RETURN m.name, r.rating, g.name", 
                  mirs = c(primate_specific_mirs$m.name, non_ps_mirs$m.name), genes = sig_tfs$tf.name)

human.count <- data.frame(dplyr::count(targets, m.name, sort = T))
head(human.count, 100)

#permutation
human.count$p.val <- NA
saveRDS(human.count, file = "./working_data/single_cell/mir_tf_permutation.rds")
} else {
  human.count <-
    readRDS(file = "./working_data/single_cell/mir_tf_permutation.rds")
}

if(ncol(human.count)<10000) {
  for (i in (ncol(human.count) + 1):10000) {
    print(i)
    perm_genes <- sample(genes, length(sig_tfs$tf.name))
    targets <- cypher(
      graph,
      "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE)
      WHERE m.name in {mirs} AND g.name in {genes} AND r.rating > 5 
      RETURN m.name, r.rating, g.name",
      mirs = c(primate_specific_mirs$m.name, non_ps_mirs$m.name), genes = perm_genes
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
  saveRDS(human.count, file = "./working_data/single_cell/mir_tf_permutation.rds")
}
tf_mirs <- human.count[human.count$p.val<.05,c("m.name", "n", "p.val")]

tf_mirs <- tf_mirs[order(tf_mirs$m.name),]
nrow(tf_mirs)
saveRDS(tf_mirs, file = "./working_data/single_cell/significant_tf_mirs.rds")

#compare with mir predictions
primate_specific_mirs <- readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <- readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")

tf_mirs[tf_mirs$m.name %in% primate_specific_mirs$m.name,]
tf_mirs[tf_mirs$m.name %in% primate_specific_mirs$m.name & tf_mirs$p.val < .001,]
nrow(tf_mirs[tf_mirs$m.name %in% primate_specific_mirs$m.name,])
tf_mirs[tf_mirs$m.name %in% non_ps_mirs$m.name & tf_mirs$p.val < .001,]
nrow(tf_mirs[tf_mirs$m.name %in% non_ps_mirs$m.name,])

#how many validated?####
val_targets <- cypher(graph, "MATCH (m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name in {mirs} AND g.name in {genes} AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.ensg", 
                  mirs = tf_mirs$m.name, genes = sig_tfs$tf.name)
dplyr::count(val_targets, g.name, sort = T)
val_targets$edge <- apply(val_targets[, c("m.name", "g.name")] , 1 , paste , collapse = "")
val_targets$pval <- tf_mirs$p.val[match(val_targets$m.name, tf_mirs$m.name)]

#mirtarbase####
mtb <- read.csv("raw_data/mirtarbase_2017_hsa.csv")
head(mtb)
g_ensg <- cypher(graph, "MATCH (g:GENE {species: 'HSA'}) WHERE g.name IN {names} RETURN g.name, g.ensg", names = mtb$Target.Gene)
g_entrez <- cypher(graph, "MATCH (g:GENE {species: 'HSA'}) WHERE g.id IN {ids} RETURN g.name, g.ensg, g.id", ids = mtb$Target.Gene..Entrez.Gene.ID.)
mtb$ensg <- g_ensg$g.ensg[match(mtb$Target.Gene, g_ensg$g.name)]
idx <- which(is.na(mtb$ensg))
mtb$ensg[idx] <- g_entrez$g.ensg[match(mtb$Target.Gene..Entrez.Gene.ID.[idx], g_entrez$g.id)]
idx <- which(is.na(mtb$ensg))
mtb$ensg[idx] <- ensg2symbol$ensembl_gene_id[match(mtb$Target.Gene[idx], ensg2symbol$hgnc_symbol)]
mtb$edge <- apply(mtb[, c("miRNA", "Target.Gene")] , 1 , paste , collapse = "")

val <- val_targets[val_targets$edge %in% mtb$edge,]
nrow(val)
nrow(val[val$pval<.001,])
