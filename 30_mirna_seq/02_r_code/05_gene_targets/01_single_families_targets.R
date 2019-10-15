#WRITE SINGLE NETWORKS####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

options(stringsAsFactors = F)

library(RNeo4j)
library(ggplot2)
library(reshape2)
library(rgexf)
library(plyr)
library(dplyr)
graph <- startGraph("http://localhost:7474/db/data/")

mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")
mir.matrix <- mir.matrix[rowMeans(mir.matrix[, 1:8]) != 0,]
pheatmap::pheatmap(asinh(as.matrix(mir.matrix[, 1:8])))
miFam <- readRDS("working_data/mir_family_table_extended.rds")

#order mir.matrix by sum of hits -> interest####
mir.matrix_bin <- apply(mir.matrix, 2, function(x) as.numeric(x != 0))
mir.matrix$de_timepoints <- rowSums(mir.matrix_bin)

#order by mean absolute countChange####
mir.matrix$de_mean_countChange <- apply(mir.matrix[, 1:8], 1, function(row) {
  mean(abs(row)[abs(row) > 0])
})

#family column####
mir.matrix$name <- rownames(mir.matrix)
mir.matrix$family <- miFam$family_name[match(mir.matrix$name, miFam$m.name)]

#network generation####
get.targeting <- function(mirs, genes = NULL) {
  if(is.null(genes)){
    targets <- cypher(graph, "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name in {mirnames}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
                    mirnames = mirs)
  } else {
    targets <- cypher(graph, "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name in {mirnames} AND g.name in {genenames}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
                      mirnames = mirs, genenames = genes)
  }
}

df2 <- readRDS("working_data/de_enriched_families.rds")

# write single or composite family networks####
write.network <- function(families) {
  message(paste("targeting for family", paste(families, collapse = " ")))
  targets <- get.targeting(unique(miFam$m.name[miFam$family_name %in% families]))
  head(targets)
  nrow(targets)
  
  de_df <- readRDS("working_data/de_counts_per_mir_unique.rds")
  
  de_df <- de_df[grep("day", de_df$time),]
  de_df <- de_df[order(abs(de_df$value), decreasing = T), ]
  de_df <- de_df[order(de_df$cell), ]
  rest <- de_df[which(duplicated(de_df$name)),]
  de_df <- de_df[-which(duplicated(de_df$name)),]
  dplyr::count(de_df, cell)
  
  targets$de <- "NONE"
  targets$de[targets$m.name %in% de_df$name[de_df$cell == "la2"]] <-
    "LA2"
  targets$de[targets$m.name %in% de_df$name[de_df$cell == "la5"]] <-
    "LA5"
  targets$de[targets$m.name %in% de_df$name[de_df$cell == "similar"]] <-
    "SIMILAR"
  targets$de[targets$m.name %in% de_df$name[de_df$cell == "opposite"]] <-
    "OPPOSITE"
  table(targets$de)
  
  nodes <-
    data.frame(id = NA, label = unique(c(targets$m.name, targets$g.name)))
  nodes$id <- seq(nrow(nodes))
  nodes$biotype <- factor(grepl("hsa", nodes$label))
  levels(nodes$biotype) <- c("gene", "miR")
  nodes$family <- NA
  idx <- which(nodes$biotype == "miR")
  
  nodes$family[idx] <- miFam$family_name[match(nodes$label[idx], miFam$m.name)]
  nodes[is.na(nodes$family) & nodes$biotype == "miR", ]
  nodes <- nodes[!(is.na(nodes$family) & nodes$biotype == "miR"), ]
  
  nodes$cell <- "NONE"
  #mirs
  nodes$cell[nodes$label %in% de_df$name[de_df$cell == "la2"]] <-
    "LA2"
  nodes$cell[nodes$label %in% de_df$name[de_df$cell == "la5"]] <-
    "LA5"
  nodes$cell[nodes$label %in% de_df$name[de_df$cell == "similar"]] <-
    "SIMILAR"
  nodes$cell[nodes$label %in% de_df$name[de_df$cell == "opposite"]] <-
    "OPPOSITE"
  #genes
  nodes$cell[nodes$label %in% targets$g.name[targets$de == "LA2"] &
               !nodes$label %in% targets$g.name[targets$de == "LA5"]] <- "LA2"
  nodes$cell[!nodes$label %in% targets$g.name[targets$de == "LA2"] &
               nodes$label %in% targets$g.name[targets$de == "LA5"]] <- "LA5"
  nodes$cell[nodes$label %in% targets$g.name[targets$de == "SIMILAR"]] <-
    "SIMILAR"
  nodes$cell[nodes$label %in% targets$g.name[targets$de == "OPPOSITE"]] <-
    "OPPOSITE"
  table(nodes$cell)
  
  nodes$bio_cell <-
    factor(apply(nodes[, c("biotype", "cell")], 1, paste, collapse = "_"))
  
  #color by bio_cell
  #female 
  fe_light <- data.frame(r = 245,g = 133,b = 255,alpha = 1)
  fe_dark <- data.frame(r = 184,g = 100,b = 191,alpha = 1)
  #male 
  ma_light <- data.frame(r = 76,g = 169,b = 255,alpha = 1)
  ma_dark <- data.frame(r = 57,g = 127,b = 191,alpha = 1)
  #similar
  si_light <- data.frame(r = 115,g = 255,b = 112,alpha = 1)
  si_dark <- data.frame(r = 86,g = 191,b = 84,alpha = 1)
  #opposite
  op_light <- data.frame(r = 255,g = 241,b = 112,alpha = 1)
  op_dark <- data.frame(r = 191,g = 181,b = 84,alpha = 1)
  #none
  no_light <- data.frame(r = 179,g = 179,b = 179,alpha = 1)
  no_dark <- data.frame(r = 134,g = 134,b = 134,alpha = 1)
  
  l <- nrow(nodes)
  color <- data.frame(r = numeric(l), g = numeric(l), b = numeric(l), alpha = numeric(l))
  for(bc in levels(nodes$bio_cell)){
    idx <- which(nodes$bio_cell == bc)
    if(bc == "miR_NONE")
      color[idx,] <- no_light
    if(bc == "gene_NONE")
      color[idx,] <- no_dark
    if(bc == "miR_LA2")
      color[idx,] <- fe_light
    if(bc == "gene_LA2")
      color[idx,] <- fe_dark
    if(bc == "miR_LA5")
      color[idx,] <- ma_light
    if(bc == "gene_LA5")
      color[idx,] <- ma_dark
    if(bc == "miR_SIMILAR")
      color[idx,] <- si_light
    if(bc == "gene_SIMILAR")
      color[idx,] <- si_dark
    if(bc == "miR_OPPOSITE")
      color[idx,] <- op_light
    if(bc == "gene_OPPOSITE")
      color[idx,] <- op_dark
  }
  
  
  nodes$countChange <-
    mir.matrix$de_mean_countChange[match(nodes$label, rownames(mir.matrix))]
  table(nodes$biotype)
  
  targets$src <- nodes$id[match(targets$m.name, nodes$label)]
  targets$tar <- nodes$id[match(targets$g.name, nodes$label)]
  edges <-
    data.frame(
      source = targets$src,
      target = targets$tar,
      Weight = targets$r.rating,
      cell = targets$de,
      NA
    )
  edges <- unique(edges)
  
  nrow(edges)
  nrow(edges[edges$Weight > 6, ])
  if(nrow(edges)>100000)
    edges <- edges[edges$Weight > 6, ]
  nrow(nodes)
  nrow(nodes[nodes$id %in% edges$source |
               nodes$id %in% edges$target, ])
  nodes <-
    nodes[nodes$id %in% edges$source | nodes$id %in% edges$target, ]
  edges <-
    edges[edges$source %in% nodes$id & edges$target %in% nodes$id, ]
  nodes <-
    nodes[nodes$id %in% edges$source | nodes$id %in% edges$target, ]
  table(nodes$biotype)
  
  
  #size by degree
  size <- numeric(nrow(nodes)) 
  names(size) <- seq(length(size))
  #mir count
  mc <- dplyr::count(edges, source, sort = T)
  idx <- which(nodes$biotype == "miR")
  size[idx] <- mc$n[match(names(size)[idx], mc$source)] / max(mc$n)
  #gene count
  gc <- dplyr::count(edges, target, sort = T)
  idx <- which(nodes$biotype == "gene")
  size[idx] <- gc$n[match(names(size)[idx], gc$target)] / max(gc$n)
  
  #mark nodes regarding single cell permutation enrichment####
  primate_specific_mirs <-
    readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
  non_ps_mirs <-
    readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
  tfs <-
    readRDS(file = "./working_data/single_cell/significant_tfs.rds")
  tf_mirs <-
    readRDS(file = "./working_data/single_cell/significant_tf_mirs.rds")
  
  nodes$ss_sig <- F
  nodes$ss_sig[nodes$label %in% c(primate_specific_mirs$m.name,
                                  non_ps_mirs$m.name,
                                  tfs$tf.name)] <- T
  nodes$tf_mir <- nodes$label %in% tf_mirs$m.name
  
  message("writing")
  fn <- paste(families, collapse = "_")
  sink(paste0("out/gephi/families/", fn, ".gexf"), type = "output")
  print(
    write.gexf(
      nodes[, c("id", "label")],
      edges[, c("source", "target")],
      edgesWeight = edges$Weight,
      nodesAtt = nodes[, 3:ncol(nodes)],
      nodesVizAtt = list(color = color, size = size^2),
      edgesAtt = edges[, 4:ncol(edges)]
    )
  )
  sink()
}

length(unique(miFam$family_name))
de_fams <- unique(miFam$family_name[miFam$de])
length(de_fams)
#single networks
lapply(de_fams, write.network)
#composites
write.network(c("mir-10", "mir-199"))
