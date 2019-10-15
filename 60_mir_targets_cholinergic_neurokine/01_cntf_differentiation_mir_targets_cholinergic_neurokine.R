#FIND DME MIR TARGETS####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")
cholinergic_genes <- readRDS("./working_data/cholinergic_genes.rds")
gene_tar <- cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")]
length(gene_tar)
#neurokines plus cholinergic genes and tfs####
tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")
nrow(tfs)

# load de mirs
mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")
mir.matrix <- mir.matrix[rowMeans(abs(mir.matrix[, 1:8]))>0,]
nrow(mir.matrix)
miFam <- readRDS("working_data/mir_family_table_extended.rds")

#order by mean absolute countChange####
mir.matrix$de_mean_countChange <- apply(mir.matrix[, 1:8], 1, function(row) {
  mean(abs(row)[abs(row) > 0])
})

# load differentially expressed genes (SCZ/BD)
load(file = "./working_data/chol_mtar_diff_genes.RData")
scz_bd_genes <- unique(chol_mtar_bd_scz$name)
bd_genes <- unique(chol_mtar_bd_scz$name[chol_mtar_bd_scz$set == "BD"])
scz_genes <- unique(chol_mtar_bd_scz$name[chol_mtar_bd_scz$set == "SCZ"])
intersect(scz_genes, bd_genes)

vis <- apply(miFam[, grep("neurokine", colnames(miFam)):ncol(miFam)], 2, as.numeric)
rownames(vis) <- miFam$m.name
vis <- vis[!duplicated(rownames(vis)),]
pheatmap::pheatmap(vis)
pheatmap::pheatmap(vis[rowSums(vis)>2,])
chol_mirs <- rownames(vis[rowSums(vis)>2,])
chol_mirs <- c(chol_mirs, miFam$m.name[miFam$family_name == "mir-132"])

# miRNet target query of differentially expressed miRs (LA2 & LA5) and genes (SCZ and BD) ####
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

targets <- get.targeting(chol_mirs, c(gene_tar, tfs$tf.name))
targets <- unique(targets)
head(targets)
hist(targets$r.rating)

ancestors <- unique(targets[,c(1,6)])
ancestors[which(duplicated(ancestors$m.name) | duplicated(ancestors$m.name, fromLast = T)),]
library(dplyr)
ancestors <- ancestors %>% group_by(m.name) %>% summarise(a.name = paste(a.name, collapse = ", "))

head(targets)
targets$g.rating <- paste0(targets$g.name, "(", targets$r.rating, ")")

#order by highest target count per miR####
#total count
nk_mir_gene <- unique(targets[, c("g.name", "r.rating", "m.name")])
nk_count <- dplyr::count(nk_mir_gene, m.name, sort = T)
head(nk_count, 20)

#perturbed count
sum(unique(targets$g.name) %in% scz_bd_genes)
nk_mir_gene <- unique(targets[targets$g.name %in% scz_bd_genes, c("g.name", "r.rating", "m.name")])
nk_count <- dplyr::count(nk_mir_gene, m.name, sort = T)
head(nk_count, 20)

nk_mir_gene <- unique(targets[,c("m.name", "g.name", "r.rating", "a.name")])

#mature or ancestor
nk_ordered <- data.frame()
for(name in nk_count$m.name){
  temp <- nk_mir_gene[nk_mir_gene$m.name == name,]
  nk_ordered <- rbind(nk_ordered, temp)
}
write.table(nk_ordered, file = "./out/neurokine_targeting_mirs_plus_chol.csv", row.names = F, sep = "\t", quote = T)

#gexf####
de_df <- readRDS("working_data/de_counts_per_mir_unique.rds")

de_df <- de_df[grep("day", de_df$time),]
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
nodes$cell[nodes$label %in% scz_genes] <- "SCZ"
nodes$cell[nodes$label %in% bd_genes] <- "BD"
nodes$cell[nodes$label %in% intersect(scz_genes, bd_genes)] <- "BOTH"
table(nodes$cell)

nodes$bio_cell <-
  factor(apply(nodes[, c("biotype", "cell")], 1, paste, collapse = "_"))

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
  if(bc == "gene_BD")
    color[idx,] <- fe_dark
  if(bc == "miR_LA5")
    color[idx,] <- ma_light
  if(bc == "gene_SCZ")
    color[idx,] <- ma_dark
  if(bc == "miR_SIMILAR")
    color[idx,] <- si_light
  if(bc == "gene_BOTH")
    color[idx,] <- si_dark
  if(bc == "miR_OPPOSITE")
    color[idx,] <- op_light
  if(bc == "gene_OPPOSITE")
    color[idx,] <- op_dark
}


#size by degree
size <- numeric(nrow(nodes)) 
names(size) <- seq(length(size))
#mir count
mc <- dplyr::count(edges, source, sort = T)
idx <- which(nodes$biotype == "miR")
# size[idx] <- mc$n[match(names(size)[idx], mc$source)] / max(mc$n)
# size[idx][is.na(size[idx])] <- min(size[idx], na.rm = T)
size[idx] <- nodes$countChange[idx] / max(nodes$countChange[idx], na.rm = T) #use countChange
size[idx][is.na(size[idx])] <- min(size[idx], na.rm = T)
hist(size[idx])
hist(sqrt(size[idx]))
size[idx] <- sqrt(size[idx])
#gene count
gc <- dplyr::count(edges, target, sort = T)
idx <- which(nodes$biotype == "gene")
size[idx] <- gc$n[match(names(size)[idx], gc$target)] / max(gc$n)
# size[idx] <- 1e-8 #ignore genes in full network

size[is.na(size)] <- min(size, na.rm = T)
hist(size)
hist(sqrt(size))


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

#mark number of targeted perturbed genes####
nodes$n_targeted_bd_scz <- nk_count$n[match(nodes$label, nk_count$m.name)]

sink(paste0("out/gephi/cholinergic/multiple_enriched.gexf"), type = "output")
print(
  write.gexf(
    nodes[, c("id", "label")],
    edges[, c("source", "target")],
    edgesWeight = edges$Weight,
    nodesAtt = nodes[, 3:ncol(nodes)],
    nodesVizAtt = list(color = color, size = size),
    edgesAtt = edges[, 4:ncol(edges)]
  )
)
sink()
