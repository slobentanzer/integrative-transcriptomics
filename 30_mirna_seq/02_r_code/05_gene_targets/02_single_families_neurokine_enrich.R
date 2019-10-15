#DE MIR FAMILIES, QUANTITATIVE####
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
miFam <- readRDS("working_data/mir_families_mirbase_v21.rds")
cholinergic_genes <- readRDS("./working_data/cholinergic_genes.rds")

#order mir.matrix by sum of hits -> interest####
mir.matrix_bin <- apply(mir.matrix, 2, function(x) as.numeric(x > 0))
mir.matrix$de_timepoints <- rowSums(mir.matrix_bin)

#order by mean absolute countChange####
mir.matrix$de_mean_countChange <- apply(mir.matrix[, 1:8], 1, function(row) {
  mean(abs(row)[abs(row) > 0])
})

#family column####
mir.matrix$family <- miFam$family_name[match(rownames(mir.matrix), miFam$m.name)]

#permutation test of targeting cholinergic/neurokine list####
#get all mirs####
allmirs <-
  cypher(graph, "MATCH (m:MIR {species: 'HSA'}) RETURN m.name")$m.name
allancs <-
  cypher(graph, "MATCH (a:ANCESTOR {species: 'HSA'}) RETURN a.name")$a.name
#>all mirs####
fams <- dplyr::count(miFam, family_name, sort = T)
fams <- fams$family_name[fams$n>4]
for(fam in fams){
  message(fam)
  targets <- cypher(graph, "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name IN {mirs} AND g.name IN {genes}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
                  mirs = miFam$m.name[miFam$family_name == fam], genes = cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")])
  targets <- unique(targets[, c("m.name", "r.rating", "g.name")])
  head(targets)
  
  #perm
  dir <- paste0("working_data/permutation/", fam)
  dir.create(dir, recursive = T)
  n <- length(list.files(dir))
  
  #count ancestors, use those to prevent similarity bias
  miFam$a.name[miFam$family_name == fam]
  m <- length(unique(miFam$a.name[miFam$family_name == fam]))
  max <- 1000
  if(n<max)
    for (i in (n + 1):max) {
      if (i %% 100 == 0)
        message(i)
      ancs <- sample(allancs, m)
      a_targets <-
        cypher(
          graph,
          "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE)
          WHERE a.name IN {ancs} AND g.name IN {genes}
          AND r.rating > 5
          RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
          ancs = ancs,
          genes = cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")]
        )
      a_targets <-
        unique(a_targets[, c("m.name", "r.rating", "g.name")])
      head(a_targets)
      count <- nrow(a_targets)
      sum <- sum(a_targets$r.rating)
      write.table(
        c(count, sum),
        file = paste0(dir, "/", i, ".txt"),
        sep = "\n",
        row.names = F,
        col.names = F
      )
    }
}
#>read####
fams0 <- vector(mode = "list", length = length(fams))
names(fams0) <- fams
fams <- fams0
for(fam in names(fams)){
  message(fam)
  #actual
  targets <- cypher(graph, "MATCH (a:ANCESTOR)-[:STEMS_FROM]-(m:MIR)-[r:ALGOSUM]-(g:GENE) 
                  WHERE m.name IN {mirs} AND g.name IN {genes}
                  AND r.rating > 5
                  RETURN m.name, r.rating, g.name, g.id, g.ensg, a.name",
                    mirs = miFam$m.name[miFam$family_name == fam], 
                    genes = cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")])
  targets <- unique(targets[, c("m.name", "r.rating", "g.name")])
  head(targets)
  if(!is.null(targets)){
    message(nrow(targets))
    dir <- paste0("working_data/permutation/", fam)
    n <- length(list.files(dir))
    null_count <- numeric()
    null_sum <- numeric()
    for (i in 1:n) {
      temp <- read.table(paste0(dir, "/", i, ".txt"), header = F)
      null_count[i] <- temp$V1[1]
      if (temp$V1[1] > 0)
        null_sum[i] <- temp$V1[2]
    }
    null_count <- sort(null_count)
    null_count <-
      cbind(null_count, data.frame(prob = 1 - abs(2 * seq(
        1:length(null_count)
      ) / length(null_count) - 1)))
    null_sum <- sort(null_sum)
    null_sum <-
      cbind(null_sum, data.frame(prob = 1 - abs(2 * seq(1:length(
        null_sum
      )) / length(null_sum) - 1)))
    
    count <- nrow(targets)
    sum <- sum(targets$r.rating)
    pval_c <-
      null_count$prob[findInterval(count, null_count$null_count)]
    pval_s <-
      null_sum$prob[findInterval(sum, null_sum$null_sum)]
    fams[[fam]] <- list(pval_c, pval_s)
  }
}
fams <- fams[!unlist(lapply(fams, is.null))]
fams[which(unlist(lapply(fams, function(x) x[[1]] < .05)))]
sig_fams <- fams[which(unlist(lapply(fams, function(x) x[[2]] < .05)))]
chol_enrich <- names(which(unlist(lapply(fams, function(x) x[[2]] < .05))))
saveRDS(sig_fams, file = "working_data/cholinergic_permutation_family_enrichment.rds")

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

# write single family networks####
gene_tar <- cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")]
for(fam in chol_enrich){
  message(fam)
  targets <- get.targeting(miFam$m.name[miFam$family_name == fam], gene_tar)
  head(targets)
  nrow(targets)
  
  de_df <- readRDS("working_data/de_counts_per_mir_unique.rds")
  
  de_df <- de_df[order(abs(de_df$value), decreasing = T), ]
  de_df <- de_df[order(de_df$cell), ]
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
    apply(nodes[, c("biotype", "cell")], 1, paste, collapse = "_")
  
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
  
  sink(paste0("out/gephi/cholinergic/", fam, ".gexf"), type = "output")
  print(
    write.gexf(
      nodes[, c("id", "label")],
      edges[, c("source", "target")],
      edgesWeight = edges$Weight,
      nodesAtt = nodes[, 3:ncol(nodes)],
      edgesAtt = edges[, 4:ncol(edges)]
    )
  )
  sink()
}

#mir-10/199 chol network####
targets <- get.targeting(miFam$m.name[miFam$family_name %in% c("mir-10", "mir-199", chol_enrich)], gene_tar)
head(targets)
nrow(targets)

de_df <- readRDS("working_data/de_counts_per_mir_unique.rds")

de_df <- de_df[order(abs(de_df$value), decreasing = T), ]
de_df <- de_df[order(de_df$cell), ]
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
  apply(nodes[, c("biotype", "cell")], 1, paste, collapse = "_")

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

sink(paste0("out/gephi/cholinergic/mir-10_mir-199_plus.gexf"), type = "output")
print(
  write.gexf(
    nodes[, c("id", "label")],
    edges[, c("source", "target")],
    edgesWeight = edges$Weight,
    nodesAtt = nodes[, 3:ncol(nodes)],
    edgesAtt = edges[, 4:ncol(edges)]
  )
)
sink()
