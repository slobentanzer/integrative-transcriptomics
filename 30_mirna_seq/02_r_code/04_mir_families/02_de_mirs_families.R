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
graph <- startGraph("http://localhost:7474/db/data/")

mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")
mir.matrix <- mir.matrix[rowMeans(mir.matrix[, 1:8]) != 0,]
pheatmap::pheatmap(asinh(as.matrix(mir.matrix[, 1:8])))
miFam <- readRDS("working_data/mir_families_mirbase_v21.rds")

#order mir.matrix by sum of hits -> interest####
mir.matrix_bin <- apply(mir.matrix, 2, function(x) as.numeric(x > 0))
mir.matrix$de_timepoints <- rowSums(mir.matrix_bin)

#order by mean absolute countChange####
mir.matrix$de_mean_countChange <- apply(mir.matrix[, 1:8], 1, function(row) {
  mean(abs(row)[abs(row) > 0])
})

#family column####
mir.matrix$name <- rownames(mir.matrix)
mir.matrix$family <- miFam$family_name[match(mir.matrix$name, miFam$m.name)]

#check top families####
#most timepoints
mir.matrix <- mir.matrix[order(mir.matrix$de_timepoints, decreasing = T),]
head(mir.matrix[, c("family", "name")], 50)
mir.matrix$tp_rank <- seq(nrow(mir.matrix))
#highest mean countchange
mir.matrix <- mir.matrix[order(mir.matrix$de_mean_countChange, decreasing = T),]
head(mir.matrix[, c("family", "name")], 50)
mir.matrix$cc_rank <- seq(nrow(mir.matrix))
#combine
mir.matrix$comb_rank <- mir.matrix$tp_rank + mir.matrix$cc_rank
mir.matrix <- mir.matrix[order(mir.matrix$comb_rank),]

#hypergeometric test in relation to family sizes####
#get all mirs####
allmirs <-
  cypher(graph, "MATCH (m:MIR {species: 'HSA'}) RETURN m.name")$m.name
# allmirs <- tolower(allmirs)
# allmirs <- gsub("-[3-5]p", "", allmirs)
# allmirs <- unique(allmirs)
#>all de mirs####
cols <- c(colnames(mir.matrix)[1:8], "family", "name")
mm_ht0 <- mir.matrix[, cols][apply(mir.matrix[, 1:8], 1, function(x) any(x != 0)),]
for(top in c(seq(10, nrow(mm_ht0), by = 10), nrow(mm_ht0))){
  mm_ht <- mm_ht0[1:top, ]
  mm_ht <- mm_ht[!duplicated(mm_ht$name), ]
  #draw 236 DE miRs from all 2588 miRs, in family vs not in family
  #per family
  # x     | family size-x
  # 236-x | 2588-family size
  
  draw <- nrow(mm_ht)
  left <- 1917 - draw
  nams <- unique(mm_ht$family)
  nams <- nams[!is.na(nams)]
  fams <- vector(mode = "list", length(nams))
  names(fams) <- sort(nams)
  for (fam in names(fams)) {
    x <- length(which(mm_ht$family == fam))
    fs <- length(which(miFam$family_name == fam))
    fams[[fam]] <-
      fisher.test(matrix(c(x, fs - x, draw - x, left + x), 2, 2), alternative = "greater")
  }
  sig_fams <- fams[unlist(lapply(fams, function(x)
    x$p.value < .05))]
  names(sig_fams)
  sig_miFam <- miFam[miFam$family_name %in% names(sig_fams), ]
  sig_miFam <- sig_miFam[sig_miFam$m.name %in% mm_ht$name,]
  if(top == 10){
    sig_fams_top10 <- sig_fams
    sig_miFam_top10 <- sig_miFam
  }
  if(top == 100){
    sig_fams_top100 <- sig_fams
    sig_miFam_top100 <- sig_miFam
  }
  if(top == nrow(mm_ht0)){
    sig_fams_all <- sig_fams
    sig_miFam_all <- sig_miFam
  }
  print(sort(table(sig_miFam$family_name)))
}
sig_fams_top10["let-7"]
sig_fams_top10["mir-10"]
sig_miFam_top10
sort(table(sig_miFam_top10$family_name))
sig_fams_top100["let-7"]
sig_fams_top100["mir-10"]
sig_miFam_top100
sort(table(sig_miFam_top100$family_name))

all_pvals <- unlist(lapply(sig_fams_all, function(x) x$p.value))
sort(all_pvals)

#>celltype specific####
#>>la2####
cols <- c(colnames(mir.matrix)[1:4], "family", "name")
mm_ht0 <- mir.matrix[, cols][apply(mir.matrix[, 1:4], 1, function(x) any(x != 0)),]
la2_nrow <- nrow(mm_ht0)
for(top in c(seq(10, nrow(mm_ht0), by = 10), nrow(mm_ht0))){
  mm_ht <- mm_ht0[1:top, ]
  mm_ht <- mm_ht[!duplicated(mm_ht$name), ]
  #draw 236 DE miRs from all 2588 miRs, in family vs not in family
  #per family
  # x     | family size-x
  # 236-x | 2588-family size
  
  draw <- nrow(mm_ht)
  left <- 1917 - draw
  nams <- unique(mm_ht$family)
  nams <- nams[!is.na(nams)]
  fams <- vector(mode = "list", length(nams))
  names(fams) <- sort(nams)
  for (fam in names(fams)) {
    x <- length(which(mm_ht$family == fam))
    fs <- length(which(miFam$family_name == fam))
    fams[[fam]] <-
      fisher.test(matrix(c(x, fs - x, draw - x, left + x), 2, 2), alternative = "greater")
  }
  sig_fams <- fams[unlist(lapply(fams, function(x)
    x$p.value < .05))]
  names(sig_fams)
  sig_miFam <- miFam[miFam$family_name %in% names(sig_fams), ]
  sig_miFam <- sig_miFam[sig_miFam$m.name %in% mm_ht$name,]
  if(top == 10){
    sig_fams_top10_la2 <- sig_fams
    sig_miFam_top10_la2 <- sig_miFam
  }
  if(top == 50){
    sig_fams_top50_la2 <- sig_fams
    sig_miFam_top50_la2 <- sig_miFam
  }
  if(top == nrow(mm_ht0)){
    sig_fams_all_la2 <- sig_fams
    sig_miFam_all_la2 <- sig_miFam
  }
  print(sort(table(sig_miFam$family_name)))
}

#>>la5####
cols <- c(colnames(mir.matrix)[5:8], "family", "name")
mm_ht0 <- mir.matrix[, cols][apply(mir.matrix[, 5:8], 1, function(x) any(x != 0)),]
nrow(mm_ht0)
for(top in c(seq(10, nrow(mm_ht0), by = 10), la2_nrow, nrow(mm_ht0))){ #la2 has 221 rows, for comparison
  mm_ht <- mm_ht0[1:top, ]
  mm_ht <- mm_ht[!duplicated(mm_ht$name), ]
  #draw 236 DE miRs from all 2588 miRs, in family vs not in family
  #per family
  # x     | family size-x
  # 236-x | 2588-family size
  
  draw <- nrow(mm_ht)
  left <- 1917 - draw
  nams <- unique(mm_ht$family)
  nams <- nams[!is.na(nams)]
  fams <- vector(mode = "list", length(nams))
  names(fams) <- sort(nams)
  for (fam in names(fams)) {
    x <- length(which(mm_ht$family == fam))
    fs <- length(which(miFam$family_name == fam))
    fams[[fam]] <-
      fisher.test(matrix(c(x, fs - x, draw - x, left + x), 2, 2), alternative = "greater")
  }
  sig_fams <- fams[unlist(lapply(fams, function(x)
    x$p.value < .05))]
  names(sig_fams)
  sig_miFam <- miFam[miFam$family_name %in% names(sig_fams), ]
  sig_miFam <- sig_miFam[sig_miFam$m.name %in% mm_ht$name,]
  if(top == 10){
    sig_fams_top10_la5 <- sig_fams
    sig_miFam_top10_la5 <- sig_miFam
  }
  if(top == 50){
    sig_fams_top50_la5 <- sig_fams
    sig_miFam_top50_la5 <- sig_miFam
  }
  if(top == la2_nrow){
    sig_fams_comp_la5 <- sig_fams
    sig_miFam_comp_la5 <- sig_miFam
  }
  if(top == nrow(mm_ht0)){
    sig_fams_all_la5 <- sig_fams
    sig_miFam_all_la5 <- sig_miFam
  }
  print(sort(table(sig_miFam$family_name)))
}

sig_fams_top10_la2["let-7"]
sig_fams_top10_la2["mir-10"]
sig_miFam_top10_la2
sig_fams_top50_la2["let-7"]
sig_fams_top50_la2["mir-10"]
sig_miFam_top50_la2
sig_fams_all_la2
sort(table(sig_miFam_all_la2$family_name))

sig_fams_top10_la5["let-7"]
sig_fams_top10_la5["mir-10"]
sig_miFam_top10_la5
sig_fams_top50_la5["let-7"]
sig_fams_top50_la5["mir-10"]
sig_miFam_top50_la5
sig_fams_all_la5
sort(table(sig_miFam_all_la5$family_name))
sig_fams_comp_la5
sort(table(sig_miFam_comp_la5$family_name))

#bar graph####
#top10 first, small break, rest
rev(sort(table(sig_miFam_top10_la2$family_name)))
rev(sort(table(sig_miFam_top10_la5$family_name)))

df <- data.frame(name = names(rev(sort(table(sig_miFam_top10_la2$family_name)))), amt = rev(sort(table(sig_miFam_top10_la2$family_name))))
df$pval_la2 <- unlist(lapply(df$name, function(x) sig_fams_top10_la2[[x]]$p.value))
df$pval_la5 <- unlist(lapply(df$name, function(x) sig_fams_top10_la5[[x]]$p.value))

#de
la2_df <- data.frame(table(sig_miFam_all_la2$family_name))
la5_df <- data.frame(table(sig_miFam_comp_la5$family_name))
#family
fam_df <- dplyr::count(miFam, family_name, sort = T)
nams <- unique(c(as.character(la2_df$Var1), as.character(la5_df$Var1)))
df2 <- data.frame(name = nams)

df2$amt_la2 <- fam_df$n[match(df2$name, fam_df$family_name)]
df2$amt_la5 <- fam_df$n[match(df2$name, fam_df$family_name)]
df2$pval_la2 <- NA
df2$pval_la5 <- NA
for(i in 1:nrow(df2)){
  nam <- df2$name[i]
  x <- sig_fams_all_la2[[nam]]$p.value
  if(!is.null(x))
    df2$pval_la2[i] <- x
  y <- sig_fams_comp_la5[[nam]]$p.value
  if(!is.null(y))
    df2$pval_la5[i] <- y
}
df2 <- df2[order(df2$pval_la2, df2$pval_la5),]
df2$pval_la2[is.na(df2$pval_la2)] <- 1
df2$log_pval_la2 <- -log10(df2$pval_la2)
df2$pval_la5[is.na(df2$pval_la5)] <- 1
df2$log_pval_la5 <- -log10(df2$pval_la5)
df2 <- df2[order(as.numeric(apply(df2[, c("log_pval_la2", "log_pval_la5")], 1, function(x) all(x != 0))), 
                 df2$amt_la2, decreasing = T),]
df2$name <- factor(df2$name, levels = df2$name)

saveRDS(df2, "working_data/de_enriched_families.rds")

df2.melt <- reshape2::melt(df2, id.vars = c("name", "amt_la2", "amt_la5"))
df2.melt$name <- factor(df2.melt$name, levels = rev(df2$name))

ggplot(df2.melt, aes(name)) + 
  geom_bar(data = subset(df2.melt, variable == "log_pval_la2"),
           aes(y = value, fill = as.factor(amt_la2)),
           stat = "identity") + 
  geom_bar(data = subset(df2.melt, variable == "log_pval_la5"),
           aes(y = -value, fill = as.factor(amt_la5)), 
           stat = "identity") + 
  coord_flip() + scale_y_reverse() +
  scale_fill_viridis_d("amt miRs") +
  xlab("miR family name") + ylab("-log10 pval")
ggsave("img/mir_families_hypergeometric_test.svg", width = 3, height = 6)

mir.matrix$la2_sig <- mir.matrix$name %in% sig_miFam_all_la2$name
mir.matrix$la5_sig <- mir.matrix$name %in% sig_miFam_all_la5$name

# #add validated to mir matrix####
# valmirs <- readRDS("working_data/validated_de_mirs_la2.rds")
# mir.matrix$validated <- rownames(mir.matrix) %in% valmirs
# 
# saveRDS(mir.matrix, file = "working_data/mir_de_matrix_countchange_extended.rds")

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

#look at significant families (all members) and mark which is DE in which cell####
nams <- miFam$m.name[miFam$family_name %in% df2$name]
targets <- get.targeting(nams)
head(targets)
# targets$short_name <- targets$m.name
# targets$short_name <- tolower(targets$short_name)
# targets$short_name <- gsub("-[3-5]p", "", targets$short_name)
# targets$short_name[grep("-[1-9]$", targets$short_name)] <- unlist(lapply(strsplit(targets$short_name[grep("-[1-9]$", targets$short_name)], "-"), function(x) paste(x[1:3], collapse = "-")))

de_df <- readRDS("working_data/de_counts_per_mir_unique.rds")

de_df <- de_df[grep("day", de_df$time),]
de_df <- de_df[order(abs(de_df$value), decreasing = T),]
de_df <- de_df[order(de_df$cell),]
de_df <- de_df[-which(duplicated(de_df$name)), ]
dplyr::count(de_df, cell)

targets$de <- "NONE"
targets$de[targets$m.name %in% de_df$name[de_df$cell == "la2"]] <- "LA2"
targets$de[targets$m.name %in% de_df$name[de_df$cell == "la5"]] <- "LA5"
targets$de[targets$m.name %in% de_df$name[de_df$cell == "similar"]] <- "SIMILAR"
targets$de[targets$m.name %in% de_df$name[de_df$cell == "opposite"]] <- "OPPOSITE"
table(targets$de)

#write complete network####
net_out <- targets
colnames(net_out) <- c("miR_name", "miRNet_score", "gene_name", "gene_entrez_id", "gene_ensg_id", "miR_precursor_name", "de_context")
write.table(net_out, row.names = F, quote = F, file = "out/DataS6.csv", sep = ";") #multiple ensgs, comma separated!####

nodes <- data.frame(id = NA, label = unique(c(targets$m.name, targets$g.name)))
nodes$id <- seq(nrow(nodes))
nodes$biotype <- factor(grepl("hsa", nodes$label))
levels(nodes$biotype) <- c("gene", "miR")
nodes$family <- NA
idx <- which(nodes$biotype == "miR")

# nams <- nodes$label[idx]
# nams <- tolower(nams)
# nams <- gsub("-[3-5]p", "", nams)
# nams[grep("-[1-9]$", nams)] <- unlist(lapply(strsplit(nams[grep("-[1-9]$", nams)], "-"), function(x) paste(x[1:3], collapse = "-")))

nodes$family[idx] <- miFam$family_name[match(nodes$label[idx], miFam$m.name)]
nodes[is.na(nodes$family) & nodes$biotype == "miR",]
nodes <- nodes[!(is.na(nodes$family) & nodes$biotype == "miR"),]

nodes$cell <- "NONE"
#mirs
nodes$cell[nodes$label %in% de_df$name[de_df$cell == "la2"]] <- "LA2"
nodes$cell[nodes$label %in% de_df$name[de_df$cell == "la5"]] <- "LA5"
nodes$cell[nodes$label %in% de_df$name[de_df$cell == "similar"]] <- "SIMILAR"
nodes$cell[nodes$label %in% de_df$name[de_df$cell == "opposite"]] <- "OPPOSITE"
#genes
nodes$cell[nodes$label %in% targets$g.name[targets$de == "LA2"] & !nodes$label %in% targets$g.name[targets$de == "LA5"]] <- "LA2"
nodes$cell[!nodes$label %in% targets$g.name[targets$de == "LA2"] & nodes$label %in% targets$g.name[targets$de == "LA5"]] <- "LA5"
nodes$cell[nodes$label %in% targets$g.name[targets$de == "SIMILAR"]] <- "SIMILAR"
nodes$cell[nodes$label %in% targets$g.name[targets$de == "OPPOSITE"]] <- "OPPOSITE"
table(nodes$cell)

nodes$bio_cell <- factor(apply(nodes[, c("biotype", "cell")], 1, paste, collapse = "_"))

nodes$countChange <- mir.matrix$de_mean_countChange[match(nodes$label, rownames(mir.matrix))]
table(nodes$biotype)

targets$src <- nodes$id[match(targets$m.name, nodes$label)]
targets$tar <- nodes$id[match(targets$g.name, nodes$label)]
edges <- data.frame(source = targets$src, target = targets$tar, Weight = targets$r.rating, cell = targets$de, NA)
edges <- unique(edges)

nrow(edges)
nrow(edges[edges$Weight > 6,])
edges <- edges[edges$Weight > 6,]
nrow(nodes)
nrow(nodes[nodes$id %in% edges$source | nodes$id %in% edges$target,])
nodes <- nodes[nodes$id %in% edges$source | nodes$id %in% edges$target,]
edges <- edges[edges$source %in% nodes$id & edges$target %in% nodes$id,]
nodes <- nodes[nodes$id %in% edges$source | nodes$id %in% edges$target,]
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
#gene count
gc <- dplyr::count(edges, target, sort = T)
idx <- which(nodes$biotype == "gene")
# size[idx] <- gc$n[match(names(size)[idx], gc$target)] / max(gc$n)
size[idx] <- 1e-8 #ignore genes in full network

size[is.na(size)] <- min(size, na.rm = T)
hist(size)
hist(sqrt(size))

#mark nodes regarding single cell permutation enrichment####
primate_specific_mirs <- readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <- readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")
tf_mirs <- readRDS(file = "./working_data/single_cell/significant_tf_mirs.rds")

nodes$ss_sig <- F
nodes$ss_sig[nodes$label %in% c(primate_specific_mirs$m.name, non_ps_mirs$m.name, tfs$tf.name)] <- T
nodes$tf_mir <- nodes$label %in% tf_mirs$m.name

sink("out/gephi/family_miRs.gexf", type = "output")
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

#summary statistics####
#mirs
#outgoing edges: number of source per node
mir_edges <- numeric()
for(mir in unique(nodes$label[nodes$biotype == "miR"])){
  mir_edges[mir] <- sum(edges$source == nodes$id[nodes$label == mir])
}
hist(mir_edges)
quantile(mir_edges)
median(mir_edges)
min(mir_edges); max(mir_edges)
mean(mir_edges); sd(mir_edges)

snodes <- nodes

miFam_stat <- miFam[miFam$m.name %in% snodes$label,]
miFam_stat <- miFam_stat[!duplicated(miFam_stat$m.name),]
for(mir in miFam_stat$m.name){
  idx <- which(miFam_stat$m.name == mir)
  miFam_stat$s_targets[idx] <- sum(edges$source[edges$cell == "SIMILAR"] %in% snodes$id[snodes$label == mir])
  miFam_stat$o_targets[idx] <- sum(edges$source[edges$cell == "OPPOSITE"] %in% snodes$id[snodes$label == mir])
  miFam_stat$f_targets[idx] <- sum(edges$source[edges$cell == "LA2"] %in% snodes$id[snodes$label == mir])
  miFam_stat$m_targets[idx] <- sum(edges$source[edges$cell == "LA5"] %in% snodes$id[snodes$label == mir])
}

miFam_stat_fam <- data.frame(miFam_stat %>% group_by(family_name) %>% dplyr::summarise(o_tars = mean(o_targets), s_tars = mean(s_targets), f_tars = mean(f_targets), m_tars = mean(m_targets)))

miFam_stat_fam$family_name <- factor(miFam_stat_fam$family_name, levels = levels(df2.melt$name))
miFam_stat_fam <- miFam_stat_fam[!is.na(miFam_stat_fam$family_name),]

mifast.melt <- melt(miFam_stat_fam)
mifast.melt$variable <- factor(mifast.melt$variable, levels = c("f_tars", "m_tars", "s_tars", "o_tars"))

ggplot(mifast.melt, aes(family_name, value, fill = variable)) + geom_bar(stat = "identity", position = "stack") +
  coord_flip()
ggsave("img/family_n_targets.svg", width = 3, height = 6)

#>statistical test####
miFam_stat$all_targets <- rowSums(miFam_stat[, c("s_targets", "o_targets", "f_targets", "m_targets")])
#sex-independent vs sex-dependent
miFam_stat$BOTH <- miFam_stat$family_name %in% df2$name[df2$pval_la2>0 & df2$pval_la5 > 0]
t.test(miFam_stat$all_targets[miFam_stat$BOTH], miFam_stat$all_targets[!miFam_stat$BOTH], alternative = "g")

#individual families against rest
tests <- list()
for(fam in unique(miFam_stat$family_name)) {
  message(fam)
  idx <- grep(paste0(fam, "$"), miFam_stat$family_name)
  if(length(idx)>1)
    tests[[fam]] <- t.test(miFam_stat$all_targets[idx], miFam_stat$all_targets[-idx], alternative = "l")
}
sig <- which(unlist(lapply(tests, function(x) x$p.value < .05)))
tests[sig]

#GO enrichment of subcluster families####
# # >let-7 5p####
# idx <- c(grep("let-7", targets$m.name), grep("miR-98", targets$m.name))
# targets_let7 <- targets[idx,]
# targets_let7 <- targets_let7[grep("5p", targets_let7$m.name),]
# targets_let7 <- targets_let7[order(targets_let7$m.name, targets_let7$g.name),]
# unique(targets_let7$m.name)
# nrow(targets_let7) == nrow(unique(targets_let7[, c("m.name", "g.name")]))
# tar_uni <- targets_let7[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_let7[idx, ]
# #different ancestors
# targets_let7 <- targets_let7[, -which(colnames(targets_let7) =="a.name")]
# targets_let7 <- unique(targets_let7)
# nrow(targets_let7) == nrow(unique(targets_let7[, c("m.name", "g.ensg")]))
# gene_counts_let7 <- dplyr::count(targets_let7, g.ensg, sort = T)
# gene_sum_let7 <- aggregate(r.rating ~ g.ensg, data = targets_let7, FUN = "sum")
# gene_sum_let7 <- gene_sum_let7[order(gene_sum_let7$r.rating, decreasing = T),]
# gene_sum_let7$count <- gene_counts_let7$n[match(gene_sum_let7$g.ensg, gene_counts_let7$g.ensg)]
# gene_sum_let7$name <- targets$g.name[match(gene_sum_let7$g.ensg, targets$g.ensg)]
# plot(gene_sum_let7$r.rating ~ gene_sum_let7$count)
# nrow(gene_sum_let7)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_let7$r.rating
# names(allGenes) <- gene_sum_let7$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# let7_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/let7_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-10 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-10"])
# targets_mir10 <- targets[idx,]
# # targets_mir10 <- targets_mir10[grep("5p", targets_mir10$m.name),]
# targets_mir10 <- targets_mir10[order(targets_mir10$m.name, targets_mir10$g.name),]
# unique(targets_mir10$m.name)
# nrow(targets_mir10) == nrow(unique(targets_mir10[, c("m.name", "g.name")]))
# tar_uni <- targets_mir10[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir10[idx, ]
# #different ancestors
# targets_mir10 <- targets_mir10[, -which(colnames(targets_mir10) =="a.name")]
# targets_mir10 <- unique(targets_mir10)
# nrow(targets_mir10) == nrow(unique(targets_mir10[, c("m.name", "g.name")]))
# gene_counts_mir10 <- dplyr::count(targets_mir10, g.ensg, sort = T)
# gene_sum_mir10 <- aggregate(r.rating ~ g.ensg, data = targets_mir10, FUN = "sum")
# gene_sum_mir10 <- gene_sum_mir10[order(gene_sum_mir10$r.rating, decreasing = T),]
# gene_sum_mir10$count <- gene_counts_mir10$n[match(gene_sum_mir10$g.ensg, gene_counts_mir10$g.ensg)]
# gene_sum_mir10$name <- targets$g.name[match(gene_sum_mir10$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir10$r.rating ~ gene_sum_mir10$count)
# nrow(gene_sum_mir10)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir10$r.rating
# names(allGenes) <- gene_sum_mir10$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir10_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir10_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-17 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-17"])
# targets_mir17 <- targets[idx,]
# # targets_mir17 <- targets_mir17[grep("5p", targets_mir17$m.name),]
# targets_mir17 <- targets_mir17[order(targets_mir17$m.name, targets_mir17$g.name),]
# unique(targets_mir17$m.name)
# nrow(targets_mir17) == nrow(unique(targets_mir17[, c("m.name", "g.name")]))
# tar_uni <- targets_mir17[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir17[idx, ]
# #different ancestors
# targets_mir17 <- targets_mir17[, -which(colnames(targets_mir17) =="a.name")]
# targets_mir17 <- unique(targets_mir17)
# nrow(targets_mir17) == nrow(unique(targets_mir17[, c("m.name", "g.name")]))
# gene_counts_mir17 <- dplyr::count(targets_mir17, g.ensg, sort = T)
# gene_sum_mir17 <- aggregate(r.rating ~ g.ensg, data = targets_mir17, FUN = "sum")
# gene_sum_mir17 <- gene_sum_mir17[order(gene_sum_mir17$r.rating, decreasing = T),]
# gene_sum_mir17$count <- gene_counts_mir17$n[match(gene_sum_mir17$g.ensg, gene_counts_mir17$g.ensg)]
# gene_sum_mir17$name <- targets$g.name[match(gene_sum_mir17$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir17$r.rating ~ gene_sum_mir17$count)
# nrow(gene_sum_mir17)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir17$r.rating
# names(allGenes) <- gene_sum_mir17$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir17_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir17_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# # >mir-199 ####
# miFam[miFam$family_name == "mir-199",]
# sort(unique(targets$a.name))
# idx <- which(targets$a.name %in% c("hsa-mir-199a-1", "hsa-mir-199a-2", "hsa-mir-199b"))
# targets_mir199 <- targets[idx,]
# # targets_mir199 <- targets_mir199[grep("5p", targets_mir199$m.name),]
# targets_mir199 <- targets_mir199[order(targets_mir199$m.name, targets_mir199$g.name),]
# unique(targets_mir199$m.name)
# nrow(targets_mir199) == nrow(unique(targets_mir199[, c("m.name", "g.name")]))
# tar_uni <- targets_mir199[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir199[idx, ]
# #different ancestors
# targets_mir199 <- targets_mir199[, -which(colnames(targets_mir199) =="a.name")]
# targets_mir199 <- unique(targets_mir199)
# nrow(targets_mir199) == nrow(unique(targets_mir199[, c("m.name", "g.name")]))
# gene_counts_mir199 <- dplyr::count(targets_mir199, g.ensg, sort = T)
# gene_sum_mir199 <- aggregate(r.rating ~ g.ensg, data = targets_mir199, FUN = "sum")
# gene_sum_mir199 <- gene_sum_mir199[order(gene_sum_mir199$r.rating, decreasing = T),]
# gene_sum_mir199$count <- gene_counts_mir199$n[match(gene_sum_mir199$g.ensg, gene_counts_mir199$g.ensg)]
# gene_sum_mir199$name <- targets$g.name[match(gene_sum_mir199$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir199$r.rating ~ gene_sum_mir199$count)
# nrow(gene_sum_mir199)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir199$r.rating
# names(allGenes) <- gene_sum_mir199$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir199_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir199_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >>mir-199 5p ####
# miFam[miFam$family_name == "mir-199",]
# sort(unique(targets$a.name))
# idx <- which(targets$a.name %in% c("hsa-mir-199a-1", "hsa-mir-199a-2", "hsa-mir-199b"))
# targets_mir199_5p <- targets[idx,]
# targets_mir199_5p <- targets_mir199_5p[grep("5p", targets_mir199_5p$m.name),]
# targets_mir199_5p <- targets_mir199_5p[order(targets_mir199_5p$m.name, targets_mir199_5p$g.name),]
# unique(targets_mir199_5p$m.name)
# nrow(targets_mir199_5p) == nrow(unique(targets_mir199_5p[, c("m.name", "g.name")]))
# tar_uni <- targets_mir199_5p[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir199_5p[idx, ]
# #different ancestors
# targets_mir199_5p <- targets_mir199_5p[, -which(colnames(targets_mir199_5p) =="a.name")]
# targets_mir199_5p <- unique(targets_mir199_5p)
# nrow(targets_mir199_5p) == nrow(unique(targets_mir199_5p[, c("m.name", "g.name")]))
# gene_counts_mir199_5p <- dplyr::count(targets_mir199_5p, g.ensg, sort = T)
# gene_sum_mir199_5p <- aggregate(r.rating ~ g.ensg, data = targets_mir199_5p, FUN = "sum")
# gene_sum_mir199_5p <- gene_sum_mir199_5p[order(gene_sum_mir199_5p$r.rating, decreasing = T),]
# gene_sum_mir199_5p$count <- gene_counts_mir199_5p$n[match(gene_sum_mir199_5p$g.ensg, gene_counts_mir199_5p$g.ensg)]
# gene_sum_mir199_5p$name <- targets$g.name[match(gene_sum_mir199_5p$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir199_5p$r.rating ~ gene_sum_mir199_5p$count)
# nrow(gene_sum_mir199_5p)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir199_5p$r.rating
# names(allGenes) <- gene_sum_mir199_5p$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir199_5p_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir199_5p_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >>mir-199 3p ####
# miFam[miFam$family_name == "mir-199",]
# sort(unique(targets$a.name))
# idx <- which(targets$a.name %in% c("hsa-mir-199a-1", "hsa-mir-199a-2", "hsa-mir-199b"))
# targets_mir199_3p <- targets[idx,]
# targets_mir199_3p <- targets_mir199_3p[grep("3p", targets_mir199_3p$m.name),]
# targets_mir199_3p <- targets_mir199_3p[order(targets_mir199_3p$m.name, targets_mir199_3p$g.name),]
# unique(targets_mir199_3p$m.name)
# nrow(targets_mir199_3p) == nrow(unique(targets_mir199_3p[, c("m.name", "g.name")]))
# tar_uni <- targets_mir199_3p[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir199_3p[idx, ]
# #different ancestors
# targets_mir199_3p <- targets_mir199_3p[, -which(colnames(targets_mir199_3p) =="a.name")]
# targets_mir199_3p <- unique(targets_mir199_3p)
# nrow(targets_mir199_3p) == nrow(unique(targets_mir199_3p[, c("m.name", "g.name")]))
# gene_counts_mir199_3p <- dplyr::count(targets_mir199_3p, g.ensg, sort = T)
# gene_sum_mir199_3p <- aggregate(r.rating ~ g.ensg, data = targets_mir199_3p, FUN = "sum")
# gene_sum_mir199_3p <- gene_sum_mir199_3p[order(gene_sum_mir199_3p$r.rating, decreasing = T),]
# gene_sum_mir199_3p$count <- gene_counts_mir199_3p$n[match(gene_sum_mir199_3p$g.ensg, gene_counts_mir199_3p$g.ensg)]
# gene_sum_mir199_3p$name <- targets$g.name[match(gene_sum_mir199_3p$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir199_3p$r.rating ~ gene_sum_mir199_3p$count)
# nrow(gene_sum_mir199_3p)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir199_3p$r.rating
# names(allGenes) <- gene_sum_mir199_3p$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir199_3p_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir199_3p_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-23 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-23"])
# targets_mir23 <- targets[idx,]
# # targets_mir23 <- targets_mir23[grep("5p", targets_mir23$m.name),]
# targets_mir23 <- targets_mir23[order(targets_mir23$m.name, targets_mir23$g.name),]
# unique(targets_mir23$m.name)
# nrow(targets_mir23) == nrow(unique(targets_mir23[, c("m.name", "g.name")]))
# tar_uni <- targets_mir23[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir23[idx, ]
# #different ancestors
# targets_mir23 <- targets_mir23[, -which(colnames(targets_mir23) =="a.name")]
# targets_mir23 <- unique(targets_mir23)
# nrow(targets_mir23) == nrow(unique(targets_mir23[, c("m.name", "g.name")]))
# gene_counts_mir23 <- dplyr::count(targets_mir23, g.ensg, sort = T)
# gene_sum_mir23 <- aggregate(r.rating ~ g.ensg, data = targets_mir23, FUN = "sum")
# gene_sum_mir23 <- gene_sum_mir23[order(gene_sum_mir23$r.rating, decreasing = T),]
# gene_sum_mir23$count <- gene_counts_mir23$n[match(gene_sum_mir23$g.ensg, gene_counts_mir23$g.ensg)]
# gene_sum_mir23$name <- targets$g.name[match(gene_sum_mir23$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir23$r.rating ~ gene_sum_mir23$count)
# nrow(gene_sum_mir23)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir23$r.rating
# names(allGenes) <- gene_sum_mir23$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir23_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir23_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# 
# # >mir-27 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-27"])
# targets_mir27 <- targets[idx,]
# # targets_mir27 <- targets_mir27[grep("5p", targets_mir27$m.name),]
# targets_mir27 <- targets_mir27[order(targets_mir27$m.name, targets_mir27$g.name),]
# unique(targets_mir27$m.name)
# nrow(targets_mir27) == nrow(unique(targets_mir27[, c("m.name", "g.name")]))
# tar_uni <- targets_mir27[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir27[idx, ]
# #different ancestors
# targets_mir27 <- targets_mir27[, -which(colnames(targets_mir27) =="a.name")]
# targets_mir27 <- unique(targets_mir27)
# nrow(targets_mir27) == nrow(unique(targets_mir27[, c("m.name", "g.name")]))
# gene_counts_mir27 <- dplyr::count(targets_mir27, g.ensg, sort = T)
# gene_sum_mir27 <- aggregate(r.rating ~ g.ensg, data = targets_mir27, FUN = "sum")
# gene_sum_mir27 <- gene_sum_mir27[order(gene_sum_mir27$r.rating, decreasing = T),]
# gene_sum_mir27$count <- gene_counts_mir27$n[match(gene_sum_mir27$g.ensg, gene_counts_mir27$g.ensg)]
# gene_sum_mir27$name <- targets$g.name[match(gene_sum_mir27$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir27$r.rating ~ gene_sum_mir27$count)
# nrow(gene_sum_mir27)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir27$r.rating
# names(allGenes) <- gene_sum_mir27$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir27_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir27_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-30 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-30"])
# targets_mir30 <- targets[idx,]
# # targets_mir30 <- targets_mir30[grep("5p", targets_mir30$m.name),]
# targets_mir30 <- targets_mir30[order(targets_mir30$m.name, targets_mir30$g.name),]
# unique(targets_mir30$m.name)
# nrow(targets_mir30) == nrow(unique(targets_mir30[, c("m.name", "g.name")]))
# tar_uni <- targets_mir30[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir30[idx, ]
# #different ancestors
# targets_mir30 <- targets_mir30[, -which(colnames(targets_mir30) =="a.name")]
# targets_mir30 <- unique(targets_mir30)
# nrow(targets_mir30) == nrow(unique(targets_mir30[, c("m.name", "g.name")]))
# gene_counts_mir30 <- dplyr::count(targets_mir30, g.ensg, sort = T)
# gene_sum_mir30 <- aggregate(r.rating ~ g.ensg, data = targets_mir30, FUN = "sum")
# gene_sum_mir30 <- gene_sum_mir30[order(gene_sum_mir30$r.rating, decreasing = T),]
# gene_sum_mir30$count <- gene_counts_mir30$n[match(gene_sum_mir30$g.ensg, gene_counts_mir30$g.ensg)]
# gene_sum_mir30$name <- targets$g.name[match(gene_sum_mir30$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir30$r.rating ~ gene_sum_mir30$count)
# nrow(gene_sum_mir30)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir30$r.rating
# names(allGenes) <- gene_sum_mir30$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir30_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir30_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# # >mir-154 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-154"])
# targets_mir154 <- targets[idx,]
# # targets_mir154 <- targets_mir154[grep("5p", targets_mir154$m.name),]
# targets_mir154 <- targets_mir154[order(targets_mir154$m.name, targets_mir154$g.name),]
# unique(targets_mir154$m.name)
# nrow(targets_mir154) == nrow(unique(targets_mir154[, c("m.name", "g.name")]))
# gene_counts_mir154 <- dplyr::count(targets_mir154, g.ensg, sort = T)
# gene_sum_mir154 <- aggregate(r.rating ~ g.ensg, data = targets_mir154, FUN = "sum")
# gene_sum_mir154 <- gene_sum_mir154[order(gene_sum_mir154$r.rating, decreasing = T),]
# gene_sum_mir154$count <- gene_counts_mir154$n[match(gene_sum_mir154$g.ensg, gene_counts_mir154$g.ensg)]
# plot(gene_sum_mir154$r.rating ~ gene_sum_mir154$count)
# nrow(gene_sum_mir154)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir154$r.rating
# names(allGenes) <- gene_sum_mir154$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir154_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir154_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-379 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-379"])
# targets_mir379 <- targets[idx,]
# # targets_mir379 <- targets_mir379[grep("5p", targets_mir379$m.name),]
# targets_mir379 <- targets_mir379[order(targets_mir379$m.name, targets_mir379$g.name),]
# unique(targets_mir379$m.name)
# nrow(targets_mir379) == nrow(unique(targets_mir379[, c("m.name", "g.name")]))
# tar_uni <- targets_mir379[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir379[idx, ]
# #different ancestors
# targets_mir379 <- targets_mir379[, -which(colnames(targets_mir379) =="a.name")]
# targets_mir379 <- unique(targets_mir379)
# nrow(targets_mir379) == nrow(unique(targets_mir379[, c("m.name", "g.name")]))
# gene_counts_mir379 <- dplyr::count(targets_mir379, g.ensg, sort = T)
# gene_sum_mir379 <- aggregate(r.rating ~ g.ensg, data = targets_mir379, FUN = "sum")
# gene_sum_mir379 <- gene_sum_mir379[order(gene_sum_mir379$r.rating, decreasing = T),]
# gene_sum_mir379$count <- gene_counts_mir379$n[match(gene_sum_mir379$g.ensg, gene_counts_mir379$g.ensg)]
# gene_sum_mir379$name <- targets$g.name[match(gene_sum_mir379$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir379$r.rating ~ gene_sum_mir379$count)
# nrow(gene_sum_mir379)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir379$r.rating
# names(allGenes) <- gene_sum_mir379$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir379_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir379_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# # >mir-500 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-500"])
# targets_mir500 <- targets[idx,]
# # targets_mir500 <- targets_mir500[grep("5p", targets_mir500$m.name),]
# targets_mir500 <- targets_mir500[order(targets_mir500$m.name, targets_mir500$g.name),]
# unique(targets_mir500$m.name)
# nrow(targets_mir500) == nrow(unique(targets_mir500[, c("m.name", "g.name")]))
# gene_counts_mir500 <- dplyr::count(targets_mir500, g.ensg, sort = T)
# gene_sum_mir500 <- aggregate(r.rating ~ g.ensg, data = targets_mir500, FUN = "sum")
# gene_sum_mir500 <- gene_sum_mir500[order(gene_sum_mir500$r.rating, decreasing = T),]
# gene_sum_mir500$count <- gene_counts_mir500$n[match(gene_sum_mir500$g.ensg, gene_counts_mir500$g.ensg)]
# plot(gene_sum_mir500$r.rating ~ gene_sum_mir500$count)
# nrow(gene_sum_mir500)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir500$r.rating
# names(allGenes) <- gene_sum_mir500$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir500_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir500_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-148 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-148"])
# targets_mir148 <- targets[idx,]
# # targets_mir148 <- targets_mir148[grep("5p", targets_mir148$m.name),]
# targets_mir148 <- targets_mir148[order(targets_mir148$m.name, targets_mir148$g.name),]
# unique(targets_mir148$m.name)
# nrow(targets_mir148) == nrow(unique(targets_mir148[, c("m.name", "g.name")]))
# tar_uni <- targets_mir148[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir148[idx, ]
# #different ancestors
# targets_mir148 <- targets_mir148[, -which(colnames(targets_mir148) =="a.name")]
# targets_mir148 <- unique(targets_mir148)
# nrow(targets_mir148) == nrow(unique(targets_mir148[, c("m.name", "g.name")]))
# gene_counts_mir148 <- dplyr::count(targets_mir148, g.ensg, sort = T)
# gene_sum_mir148 <- aggregate(r.rating ~ g.ensg, data = targets_mir148, FUN = "sum")
# gene_sum_mir148 <- gene_sum_mir148[order(gene_sum_mir148$r.rating, decreasing = T),]
# gene_sum_mir148$count <- gene_counts_mir148$n[match(gene_sum_mir148$g.ensg, gene_counts_mir148$g.ensg)]
# gene_sum_mir148$name <- targets$g.name[match(gene_sum_mir148$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir148$r.rating ~ gene_sum_mir148$count)
# nrow(gene_sum_mir148)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir148$r.rating
# names(allGenes) <- gene_sum_mir148$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir148_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir148_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# 
# # >mir-329 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-329"])
# targets_mir329 <- targets[idx,]
# # targets_mir329 <- targets_mir329[grep("5p", targets_mir329$m.name),]
# targets_mir329 <- targets_mir329[order(targets_mir329$m.name, targets_mir329$g.name),]
# unique(targets_mir329$m.name)
# nrow(targets_mir329) == nrow(unique(targets_mir329[, c("m.name", "g.name")]))
# tar_uni <- targets_mir329[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir329[idx, ]
# #different ancestors
# targets_mir329 <- targets_mir329[, -which(colnames(targets_mir329) =="a.name")]
# targets_mir329 <- unique(targets_mir329)
# nrow(targets_mir329) == nrow(unique(targets_mir329[, c("m.name", "g.name")]))
# gene_counts_mir329 <- dplyr::count(targets_mir329, g.ensg, sort = T)
# gene_sum_mir329 <- aggregate(r.rating ~ g.ensg, data = targets_mir329, FUN = "sum")
# gene_sum_mir329 <- gene_sum_mir329[order(gene_sum_mir329$r.rating, decreasing = T),]
# gene_sum_mir329$count <- gene_counts_mir329$n[match(gene_sum_mir329$g.ensg, gene_counts_mir329$g.ensg)]
# gene_sum_mir329$name <- targets$g.name[match(gene_sum_mir329$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir329$r.rating ~ gene_sum_mir329$count)
# nrow(gene_sum_mir329)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir329$r.rating
# names(allGenes) <- gene_sum_mir329$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir329_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir329_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# # >mir-374 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-374"])
# targets_mir374 <- targets[idx,]
# # targets_mir374 <- targets_mir374[grep("5p", targets_mir374$m.name),]
# targets_mir374 <- targets_mir374[order(targets_mir374$m.name, targets_mir374$g.name),]
# unique(targets_mir374$m.name)
# nrow(targets_mir374) == nrow(unique(targets_mir374[, c("m.name", "g.name")]))
# tar_uni <- targets_mir374[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir374[idx, ]
# #different ancestors
# targets_mir374 <- targets_mir374[, -which(colnames(targets_mir374) =="a.name")]
# targets_mir374 <- unique(targets_mir374)
# nrow(targets_mir374) == nrow(unique(targets_mir374[, c("m.name", "g.name")]))
# gene_counts_mir374 <- dplyr::count(targets_mir374, g.ensg, sort = T)
# gene_sum_mir374 <- aggregate(r.rating ~ g.ensg, data = targets_mir374, FUN = "sum")
# gene_sum_mir374 <- gene_sum_mir374[order(gene_sum_mir374$r.rating, decreasing = T),]
# gene_sum_mir374$count <- gene_counts_mir374$n[match(gene_sum_mir374$g.ensg, gene_counts_mir374$g.ensg)]
# gene_sum_mir374$name <- targets$g.name[match(gene_sum_mir374$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir374$r.rating ~ gene_sum_mir374$count)
# nrow(gene_sum_mir374)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir374$r.rating
# names(allGenes) <- gene_sum_mir374$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir374_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir374_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)
# 
# # >mir-515 ####
# idx <- which(targets$a.name %in% miFam$name[miFam$family_name == "mir-515"])
# targets_mir515 <- targets[idx,]
# # targets_mir515 <- targets_mir515[grep("5p", targets_mir515$m.name),]
# targets_mir515 <- targets_mir515[order(targets_mir515$m.name, targets_mir515$g.name),]
# unique(targets_mir515$m.name)
# nrow(targets_mir515) == nrow(unique(targets_mir515[, c("m.name", "g.name")]))
# tar_uni <- targets_mir515[, c("m.name", "g.name")]
# idx <- duplicated(apply(tar_uni, 1, paste, collapse = "")) | duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
# targets_mir515[idx, ]
# #different ancestors
# targets_mir515 <- targets_mir515[, -which(colnames(targets_mir515) =="a.name")]
# targets_mir515 <- unique(targets_mir515)
# nrow(targets_mir515) == nrow(unique(targets_mir515[, c("m.name", "g.name")]))
# gene_counts_mir515 <- dplyr::count(targets_mir515, g.ensg, sort = T)
# gene_sum_mir515 <- aggregate(r.rating ~ g.ensg, data = targets_mir515, FUN = "sum")
# gene_sum_mir515 <- gene_sum_mir515[order(gene_sum_mir515$r.rating, decreasing = T),]
# gene_sum_mir515$count <- gene_counts_mir515$n[match(gene_sum_mir515$g.ensg, gene_counts_mir515$g.ensg)]
# gene_sum_mir515$name <- targets$g.name[match(gene_sum_mir515$g.ensg, targets$g.ensg)]
# plot(gene_sum_mir515$r.rating ~ gene_sum_mir515$count)
# nrow(gene_sum_mir515)
# 
# library(topGO)
# geneID2GO <- readRDS(file = "./working_data/geneID2GO_all.rds")
# 
# # BPTerms <- ls(GOBPTerm)
# # head(BPTerms)
# #ontology
# ontology <- "BP"
# #description
# description <- "miR targets in LA cells"
# #allGenes list
# allGenes <- gene_sum_mir515$r.rating
# names(allGenes) <- gene_sum_mir515$g.ensg
# #geneSelection function
# hist(allGenes)
# quantile(allGenes, seq(0,1,.01))
# x <- quantile(allGenes, seq(0,1,.01))[91]
# topGenes <- function(allScore) {
#   return(allScore >= x)
# }
# #prune GO terms
# nodeSize <- 10
# 
# GOdata <- new("topGOdata",
#               description = description,
#               ontology = ontology,
#               allGenes = allGenes,
#               geneSel = topGenes,
#               annot = annFUN.gene2GO,
#               nodeSize = nodeSize,
#               gene2GO = geneID2GO)
# 
# ## >fisher test ###
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# resultFisher
# 
# # >weight test ###
# test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
# resultWeight <- getSigGroups(GOdata, test.stat)
# resultWeight
# 
# ## >total ###
# allRes <- GenTable(GOdata, classic = resultFisher, 
#                    # KS = resultKS, 
#                    weight = resultWeight, 
#                    # elim = resultElim,
#                    orderBy = "weight", ranksOf = "classic", topNodes = 100)
# # allRes <- allRes[order(allRes$Annotated, decreasing = T),]
# mir515_sig <- allRes[allRes$weight<.05,]
# 
# write.table(allRes[allRes$weight<.05,], file = "./out/mir515_mir_targets_GO_enrich.csv",
#             quote = T, row.names = F)


# >GENERIC ####
source("GO.R")
get.GO.family <- function(family_name, suffix = NULL, limit = 3000){
  idx <-
    which(targets$m.name %in% miFam$m.name[miFam$family_name == family_name])
  targets_mirfam <- targets[idx, ]
  if(!is.null(suffix)) {
    targets_mirfam <- targets_mirfam[grep(suffix, targets_mirfam$m.name),]
  }
  targets_mirfam <-
    targets_mirfam[order(targets_mirfam$m.name, targets_mirfam$g.name), ]
  unique(targets_mirfam$m.name)
  nrow(targets_mirfam) == nrow(unique(targets_mirfam[, c("m.name", "g.name")]))
  # tar_uni <- targets_mirfam[, c("m.name", "g.name")]
  # idx <-
  #   duplicated(apply(tar_uni, 1, paste, collapse = "")) |
  #   duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
  # targets_mirfam[idx,]
  # #different ancestors
  # targets_mirfam <-
  #   targets_mirfam[,-which(colnames(targets_mirfam) == "a.name")]
  targets_mirfam <- unique(targets_mirfam)
  nrow(targets_mirfam) == nrow(unique(targets_mirfam[, c("m.name", "g.name")]))
  gene_counts_mirfam <- dplyr::count(targets_mirfam, g.ensg, sort = T)
  gene_sum_mirfam <-
    aggregate(r.rating ~ g.ensg, data = targets_mirfam, FUN = "sum")
  gene_sum_mirfam <-
    gene_sum_mirfam[order(gene_sum_mirfam$r.rating, decreasing = T), ]
  gene_sum_mirfam$count <-
    gene_counts_mirfam$n[match(gene_sum_mirfam$g.ensg, gene_counts_mirfam$g.ensg)]
  gene_sum_mirfam$name <-
    targets$g.name[match(gene_sum_mirfam$g.ensg, targets$g.ensg)]
  plot(gene_sum_mirfam$r.rating ~ gene_sum_mirfam$count)
  nrow(gene_sum_mirfam)
  
  scores <- gene_sum_mirfam$r.rating
  ensgs <- gene_sum_mirfam$g.ensg
  if(length(scores)>limit) {
    scores <- scores[1:limit]
    ensgs <- ensgs[1:limit]
  }
  
  res <- get.GO(scores, ensgs)
  
  write.table(res[[1]],
              file = paste0("./out/", family_name, "_mir_targets_GO_enrich.csv"),
              quote = T,
              row.names = F)
  
  return(list(mir_names = unique(targets_mirfam$m.name), go_table = res[[1]]))
}

mir_fam_go_return <- lapply(df2$name, get.GO.family)
mir_fam_go_return[["let-7"]] <- get.GO.family("let-7", suffix = "5p")

names(mir_fam_go_return) <- df2$name
lapply(mir_fam_go_return, "[[", 1)

#>save####
mir_fam_go_list <- lapply(mir_fam_go_return, "[[", 2)

saveRDS(mir_fam_go_list,
     file = "working_data/mir_family_GO_enrich_list.rds")

# #GO enrichment of subcluster families - only DE members####

# >GENERIC ####
get.GO.families.DE <- function(family_name, suffix = NULL, limit = 3000){
  idx <-
    which(targets$m.name %in% miFam$m.name[miFam$family_name == family_name])
  targets_mirfam <- targets[idx, ]
  targets_mirfam <- targets_mirfam[targets_mirfam$m.name %in% rownames(mir.matrix),]
  if(!is.null(suffix)) {
    targets_mirfam <- targets_mirfam[grep(suffix, targets_mirfam$m.name),]
  }
  targets_mirfam <-
    targets_mirfam[order(targets_mirfam$m.name, targets_mirfam$g.name), ]
  unique(targets_mirfam$m.name)
  nrow(targets_mirfam) == nrow(unique(targets_mirfam[, c("m.name", "g.name")]))
  tar_uni <- targets_mirfam[, c("m.name", "g.name")]
  idx <-
    duplicated(apply(tar_uni, 1, paste, collapse = "")) |
    duplicated(apply(tar_uni, 1, paste, collapse = ""), fromLast = T)
  targets_mirfam[idx,]
  #different ancestors
  targets_mirfam <-
    targets_mirfam[,-which(colnames(targets_mirfam) == "a.name")]
  targets_mirfam <- unique(targets_mirfam)
  nrow(targets_mirfam) == nrow(unique(targets_mirfam[, c("m.name", "g.name")]))
  gene_counts_mirfam <- dplyr::count(targets_mirfam, g.ensg, sort = T)
  gene_sum_mirfam <-
    aggregate(r.rating ~ g.ensg, data = targets_mirfam, FUN = "sum")
  gene_sum_mirfam <-
    gene_sum_mirfam[order(gene_sum_mirfam$r.rating, decreasing = T), ]
  gene_sum_mirfam$count <-
    gene_counts_mirfam$n[match(gene_sum_mirfam$g.ensg, gene_counts_mirfam$g.ensg)]
  gene_sum_mirfam$name <-
    targets$g.name[match(gene_sum_mirfam$g.ensg, targets$g.ensg)]
  plot(gene_sum_mirfam$r.rating ~ gene_sum_mirfam$count)
  nrow(gene_sum_mirfam)
  
  scores <- gene_sum_mirfam$r.rating
  ensgs <- gene_sum_mirfam$g.ensg
  if(length(scores)>limit) {
    scores <- scores[1:limit]
    ensgs <- ensgs[1:limit]
  }
  
  res <- get.GO(scores, ensgs)
  
  write.table(res[[1]],
              file = paste0("./out/", family_name, "_DE_mir_targets_GO_enrich.csv"),
              quote = T,
              row.names = F)
  
  return(list(mir_names = unique(targets_mirfam$m.name), go_table = res[[1]]))
}

mir_fam_go_return_de <- lapply(df2$name, get.GO.families.DE)
names(mir_fam_go_return_de) <- df2$name
lapply(mir_fam_go_return_de, "[[", 1)


#>save####
mir_fam_go_list_de <- lapply(mir_fam_go_return_de, "[[", 2)

saveRDS(mir_fam_go_list_de,
        file = "working_data/mir_family_GO_enrich_list_de.rds")

#intersect of terms####
mir_fam_go_list <- readRDS(file = "working_data/mir_family_GO_enrich_list.rds")
for(nam in names(mir_fam_go_list)){
  message(nam)
  print(mir_fam_go_list[[nam]])
  readline("ENTER")
}


library(venn)
ids <- sort(unique(plyr::ldply(mir_fam_go_list)$GO.ID))
id_matrix <- matrix(data = NA, nrow = length(ids), ncol = length(mir_fam_go_list))
colnames(id_matrix) <- names(mir_fam_go_list)
rownames(id_matrix) <- ids
for(id in ids){
  for(nam in colnames(id_matrix))
    id_matrix[id, nam] <- as.numeric(id %in% mir_fam_go_list[[nam]]$GO.ID)
}
venn(data.frame(id_matrix[, 1:4]))
nrow(id_matrix[rowSums(id_matrix)>1,])
nrow(id_matrix[rowSums(id_matrix)==1,])

#what is overlap?
overlap <- rownames(id_matrix)[rowSums(id_matrix)>1]
ol_frame <- data.frame()
excl_mir_fam_go_list <- mir_fam_go_list
for(nam in names(mir_fam_go_list)){
  temp <- mir_fam_go_list[[nam]]
  ol_frame <- rbind(ol_frame, temp[temp$GO.ID %in% overlap,])
  excl_mir_fam_go_list[[nam]] <- temp[!temp$GO.ID %in% overlap,]
}
ol_frame <- unique(ol_frame)
ol_frame$n <- unlist(lapply(ol_frame$GO.ID, function(id) sum(as.numeric(unlist(lapply(mir_fam_go_list, function(df) id %in% df$GO.ID))))))
ol_frame <- ol_frame[order(ol_frame$n, ol_frame$GO.ID, decreasing = T),]
ol_frame[!duplicated(ol_frame$GO.ID),]

#compare to figure 1 results####
#master
res.scz.bd <- read.table("out/sheet1.csv", header = T)
#non-sex
res.scz.minus.bd <- read.table("out/sheet2.csv", header = T)
res.bd.minus.scz <- read.table("out/sheet3.csv", header = T)
res.bd.sex <- read.table("out/sheet8.csv", header = T)
res.scz.sex <- read.table("out/sheet11.csv", header = T)

#disease differences per sex A + B
#A
#right
res.scz.minus.bd_m <- read.table("out/sheet4.csv", header = T)
#left
res.scz.minus.bd_f <- read.table("out/sheet5.csv", header = T)
#B
#right
res.bd.minus.scz_m <- read.table("out/sheet6.csv", header = T)
#left
res.bd.minus.scz_f <- read.table("out/sheet7.csv", header = T)

#sex differences per disease C + D
#C
#left
res.bd.m.minus.f <- read.table("out/sheet10.csv", header = T)
#right
res.scz.m.minus.f <- read.table("out/sheet13.csv", header = T)
#D
#left
res.bd.f.minus.m <- read.table("out/sheet9.csv", header = T)
#right
res.scz.f.minus.m <- read.table("out/sheet12.csv", header = T)

meta_go_list <- list(res.scz.minus.bd_f, res.scz.minus.bd_m, res.bd.minus.scz_f, res.bd.minus.scz_m,
                     res.bd.m.minus.f, res.scz.m.minus.f, res.bd.f.minus.m, res.scz.f.minus.m)
names(meta_go_list) <- c("AL", "AR", "BL", "BR", "CL", "CR", "DL", "DR")

mir_fam_go_list <- lapply(mir_fam_go_list, function(x) {
  temp <- x
  colnames(temp) <- gsub(" ", ".", colnames(temp))
  return(temp)
})
excl_mir_fam_go_list <- lapply(excl_mir_fam_go_list, function(x) {
  temp <- x
  colnames(temp) <- gsub(" ", ".", colnames(temp))
  return(temp)
})

comb_list <- c(mir_fam_go_list, meta_go_list)

#curate go_table for interesting terms####
go_table <- unique(plyr::ldply(comb_list)[, c("GO.ID", "Term")])
go_table <- go_table[order(go_table$GO.ID),]
rownames(go_table) <- seq(nrow(go_table))
head(go_table)
{
go_table <- go_table[!grepl("autophagosome", go_table$Term),]
go_table <- go_table[!grepl("cyclin-dependent", go_table$Term),]
go_table <- go_table[!grepl("mitotic", go_table$Term),]
go_table <- go_table[!grepl("MAPK", go_table$Term),]
go_table <- go_table[!grepl("polyubiq", go_table$Term),]
go_table <- go_table[!grepl("lung", go_table$Term),]
go_table <- go_table[!grepl("polyubiq", go_table$Term),]
go_table <- go_table[!grepl("DNA damage", go_table$Term),]
go_table <- go_table[!grepl("ossification", go_table$Term),]
go_table <- go_table[!grepl("minerali", go_table$Term),]
go_table <- go_table[!grepl("respiratory", go_table$Term),]
go_table <- go_table[!grepl("cholesterol", go_table$Term),]
go_table <- go_table[!grepl("megakaryo", go_table$Term),]
go_table <- go_table[!grepl("erythrocyte", go_table$Term),]
go_table <- go_table[!grepl("platelet", go_table$Term),]
go_table <- go_table[!grepl("sphingolip", go_table$Term),]
go_table <- go_table[!grepl("endocytosis", go_table$Term),]
go_table <- go_table[!grepl("cell polarity", go_table$Term),]
go_table <- go_table[!grepl("second-messenger", go_table$Term),]
go_table <- go_table[!grepl("pancreas", go_table$Term),]
go_table <- go_table[!grepl("actin", go_table$Term),]
go_table <- go_table[!grepl("actomyosin", go_table$Term),]
go_table <- go_table[!grepl("endoplasmic", go_table$Term),]
go_table <- go_table[!grepl("epithelial", go_table$Term),]
go_table <- go_table[!grepl("regulation of transcription", go_table$Term),]
go_table <- go_table[!grepl("microtubule", go_table$Term),]
go_table <- go_table[!grepl("nuclear division", go_table$Term),]
go_table <- go_table[!grepl("nuclear-transcribed", go_table$Term),]
go_table <- go_table[!grepl("reactive oxygen", go_table$Term),]
go_table <- go_table[!grepl("non-recombinational", go_table$Term),]
go_table <- go_table[!grepl("retinoid", go_table$Term),]
go_table <- go_table[!grepl("osteoblast", go_table$Term),]
go_table <- go_table[!grepl("eye", go_table$Term),]
go_table <- go_table[!grepl("metanephros", go_table$Term),]
go_table <- go_table[!grepl("ameboidal", go_table$Term),]
go_table <- go_table[!grepl("glucose", go_table$Term),]
go_table <- go_table[!grepl("embryonic", go_table$Term),]
go_table <- go_table[!grepl("primary germ layer", go_table$Term),]
go_table <- go_table[!grepl("mesoderm", go_table$Term),]
go_table <- go_table[!grepl("cell fate", go_table$Term),]
go_table <- go_table[!grepl("kidney development", go_table$Term),]
go_table <- go_table[!grepl("epithelial", go_table$Term),]
go_table <- go_table[!grepl("protein phosphory", go_table$Term),]
go_table <- go_table[!grepl("hair follicle", go_table$Term),]
go_table <- go_table[!grepl("endothelial", go_table$Term),]
go_table <- go_table[!grepl("liver development", go_table$Term),]
go_table <- go_table[!grepl("protein insertion into mitochondrial", go_table$Term),]
go_table <- go_table[!grepl("heart looping", go_table$Term),]
go_table <- go_table[!grepl("regulation of cell-matrix", go_table$Term),]
go_table <- go_table[!grepl("epithelium", go_table$Term),]
go_table <- go_table[!grepl("chondrocyte", go_table$Term),]
go_table <- go_table[!grepl("regionalization", go_table$Term),]
go_table <- go_table[!grepl("heart morphogenesis", go_table$Term),]
go_table <- go_table[!grepl("systemic arterial blood", go_table$Term),]
go_table <- go_table[!grepl("outflow tract", go_table$Term),]
go_table <- go_table[!grepl("endothelium", go_table$Term),]
go_table <- go_table[!grepl("heart", go_table$Term),]
go_table <- go_table[!grepl("ventricular", go_table$Term),]
go_table <- go_table[!grepl("endocardial", go_table$Term),]
go_table <- go_table[!grepl("cardiac", go_table$Term),]
go_table <- go_table[!grepl("axis elongation", go_table$Term),]
go_table <- go_table[!grepl("glycosaminoglycan", go_table$Term),]
go_table <- go_table[!grepl("aldehyde", go_table$Term),]
go_table <- go_table[!grepl("glycolytic", go_table$Term),]
go_table <- go_table[!grepl("DNA", go_table$Term),]
go_table <- go_table[!grepl("chromatin", go_table$Term),]
go_table <- go_table[!grepl("protein folding", go_table$Term),]
go_table <- go_table[!grepl("protein kinase", go_table$Term),]
go_table <- go_table[!grepl("deacetylation", go_table$Term),]
go_table <- go_table[!grepl("glycosylation", go_table$Term),]
go_table <- go_table[!grepl("ubiquitin", go_table$Term),]
go_table <- go_table[!grepl("protein targeting", go_table$Term),]
go_table <- go_table[!grepl("fatty acid", go_table$Term),]
go_table <- go_table[!grepl("triglyceride", go_table$Term),]
go_table <- go_table[!grepl("phospholipid", go_table$Term),]
go_table <- go_table[!grepl("phosphatidylinositol", go_table$Term),]
go_table <- go_table[!grepl("ceramide", go_table$Term),]
go_table <- go_table[!grepl("terpenoid", go_table$Term),]
go_table <- go_table[!grepl("nitric oxide", go_table$Term),]
go_table <- go_table[!grepl("calcium ion", go_table$Term),]
go_table <- go_table[!grepl("amino acid transport", go_table$Term),]
go_table <- go_table[!grepl("Golgi", go_table$Term),]
go_table <- go_table[!grepl("nucleocytoplasmic", go_table$Term),]
go_table <- go_table[!grepl("muscle contraction", go_table$Term),]
go_table <- go_table[!grepl("nucleus organization", go_table$Term),]
go_table <- go_table[!grepl("cytoskeleton", go_table$Term),]
go_table <- go_table[!grepl("endosome", go_table$Term),]
go_table <- go_table[!grepl("vacuolar", go_table$Term),]
go_table <- go_table[!grepl("meiosis", go_table$Term),]
go_table <- go_table[!grepl("body fluid secretion", go_table$Term),]
go_table <- go_table[!grepl("coagulation", go_table$Term),]
go_table <- go_table[!grepl("perception", go_table$Term),]
go_table <- go_table[!grepl("blood circulation", go_table$Term),]
go_table <- go_table[!grepl("regulation of cell", go_table$Term),]
go_table <- go_table[!grepl("axon ensheathment", go_table$Term),]
go_table <- go_table[!grepl("epidermis", go_table$Term),]
go_table <- go_table[!grepl("glycoprotein", go_table$Term),]
go_table <- go_table[!grepl("animal organ", go_table$Term),]
go_table <- go_table[!grepl("catabolic process", go_table$Term),]
go_table <- go_table[!grepl("anterior/posterior", go_table$Term),]
go_table <- go_table[!grepl("radiation", go_table$Term),]
go_table <- go_table[!grepl("endomembrane", go_table$Term),]
go_table <- go_table[!grepl("macromolecule", go_table$Term),]
go_table <- go_table[!grepl("autophagy", go_table$Term),]
go_table <- go_table[!grepl("macromolecule", go_table$Term),]
go_table <- go_table[!grepl("organelle", go_table$Term),]
go_table <- go_table[!grepl("fibroblast", go_table$Term),]
go_table <- go_table[!grepl("ammonium", go_table$Term),]
go_table <- go_table[!grepl("bicarbonate", go_table$Term),]
go_table <- go_table[!grepl("monocarboxylic acid", go_table$Term),]
go_table <- go_table[!grepl("peptidyl", go_table$Term),]
go_table <- go_table[!grepl("phenol-containing", go_table$Term),]
go_table <- go_table[!grepl("cardiocyte", go_table$Term),]
go_table <- go_table[!grepl("pigmentation", go_table$Term),]
go_table <- go_table[!grepl("carbohydrate", go_table$Term),]
go_table <- go_table[!grepl("ear", go_table$Term),]
go_table <- go_table[!grepl("skin", go_table$Term),]
go_table <- go_table[!grepl("donto", go_table$Term),]
go_table <- go_table[!grepl("fat cell", go_table$Term),]
go_table <- go_table[!grepl("osteoclast", go_table$Term),]
go_table <- go_table[!grepl("keratinocyte", go_table$Term),]
go_table <- go_table[!grepl("bone remodeling", go_table$Term),]
go_table <- go_table[!grepl("lipid biosynthetic", go_table$Term),]
go_table <- go_table[!grepl("iron ion", go_table$Term),]
go_table <- go_table[!grepl("adipose tissue", go_table$Term),]
go_table <- go_table[!grepl("somite", go_table$Term),]
go_table <- go_table[!grepl("cochlea", go_table$Term),]
go_table <- go_table[!grepl("COPII", go_table$Term),]
go_table <- go_table[!grepl("ciliary basal", go_table$Term),]
go_table <- go_table[!grepl("import into", go_table$Term),]
go_table <- go_table[!grepl("response to ketone", go_table$Term),]
go_table <- go_table[!grepl("semi-lunar", go_table$Term),]
go_table <- go_table[!grepl("roof of mouth", go_table$Term),]
go_table <- go_table[!grepl("muscle", go_table$Term),]
go_table <- go_table[!grepl("anatomical", go_table$Term),]
go_table <- go_table[!grepl("hexose", go_table$Term),]
go_table <- go_table[!grepl("myoblast", go_table$Term),]
go_table <- go_table[!grepl("insulin", go_table$Term),]
go_table <- go_table[!grepl("XXX", go_table$Term),]
}
pertinent_ids <- go_table$GO.ID

#print overlap####
go_table <- unique(plyr::ldply(comb_list)[, c("GO.ID", "Term")])

ids <- sort(unique(plyr::ldply(comb_list)$GO.ID))
id_matrix <- matrix(data = NA, nrow = length(ids), ncol = length(comb_list))
colnames(id_matrix) <- names(comb_list)
rownames(id_matrix) <- ids
for(id in ids){
  for(nam in colnames(id_matrix))
    id_matrix[id, nam] <- as.numeric(id %in% comb_list[[nam]]$GO.ID)
}

combs <- combn(seq(length(comb_list)), 2)
n <- ncol(combs)
overlap_list <- vector("list", n)
for(i in 1:n){
  cols <- combs[, i]
  # venn(data.frame(id_matrix[, cols]))
  overlap_list[[i]] <- go_table[go_table$GO.ID %in% rownames(id_matrix[rowSums(id_matrix[, cols]) > 1, cols]),]
}
names(overlap_list) <- apply(combs, 2, function(x) paste(colnames(id_matrix)[x[1]], colnames(id_matrix)[x[2]], sep = "_"))
overlap_list <- overlap_list[unlist(lapply(overlap_list, function(x) nrow(x)>0))]

overlap_frame <- plyr::ldply(overlap_list)
overlap_count <- vector("numeric")
for(nam in names(mir_fam_go_list)){
  message(nam)
  temp <- overlap_frame[grep(paste0(nam, "_"), overlap_frame$.id),]
  # temp <- temp[temp$GO.ID %in% pertinent_ids,]
  overlap_count[nam] <- nrow(temp)
  print(unique(temp[grep("[A-D]", temp$.id),]))
  # readline("ENTER")
}
overlap_count <- data.frame(name = names(overlap_count), count = overlap_count)
rownames(overlap_count) <- seq(nrow(overlap_count))
overlap_count$enriched <- c(rep("BOTH", 5), rep("ONE", 12))
mean(overlap_count$count[overlap_count$enriched == "BOTH"])
mean(overlap_count$count[overlap_count$enriched == "ONE"])
t.test(x = overlap_count$count[overlap_count$enriched == "BOTH"],
       y = overlap_count$count[overlap_count$enriched == "ONE"],
       alternative = "g")

#print pertinent terms####
length(unique(c(plyr::ldply(mir_fam_go_list)$GO.ID, plyr::ldply(mir_fam_go_list_de)$GO.ID)))


for(nam in names(mir_fam_go_list)){
  temp <- mir_fam_go_list[[nam]]
  print(temp[temp$GO.ID %in% pertinent_ids,])
  message(nam) 
  readline("ENTER")
}

#curate topics####
rows.list <- list()
rows.list_de <- list()
#topic names: 
topics <- c("sex", "inflammation", "t_cells", "b_cells", "circadian", "neuron_death", "aging", "cytokines", 
            "neurotrophins", "splicing", "microRNA", "Z new topic?", "leukocytes", "nfkb", "substance_response", "neurodevelopment", "plasticity", "apoptosis")
for(nam in names(mir_fam_go_list_de)){
  #for each family
  message(nam)
  start <- readline("Start?")
  if(start == "y") {
    count <- 1
    while (count <= length(topics)) {
      n_topic <- NULL
      topics <- sort(topics)
      topic <- topics[count]
      message(topic)
      #for each topic
      if (topic == "Z new topic?") {
        n_topic <- readline("ENTER TOPIC: ")
        if (nchar(n_topic) > 0)
          topics <- c(n_topic, topics)
      }
      
      if (is.null(n_topic)) {
        if (is.null(rows.list[[topic]]))
          rows.list[[topic]] <- list()
        if (is.null(rows.list_de[[topic]]))
          rows.list_de[[topic]] <- list()
        
        #print and enter
        temp <- mir_fam_go_list[[nam]]
        print(temp[temp$GO.ID %in% pertinent_ids,])
        message(nam)
        nums <- readline("ENTER ROWS: ")
        if (nchar(nums) > 0)
          rows.list[[topic]][[nam]] <-
          as.numeric(unlist(strsplit(nums, ",")))
        
        temp <- mir_fam_go_list_de[[nam]]
        print(temp[temp$GO.ID %in% pertinent_ids, ])
        message(nam)
        numsDE <- readline("ENTER ROWS (DE): ")
        if (nchar(numsDE) > 0)
          rows.list_de[[topic]][[nam]] <-
          as.numeric(unlist(strsplit(numsDE, ",")))
        count <- count + 1
      } else {
        if (nchar(n_topic) > 0) {
          topic <- n_topic
          if (is.null(rows.list[[topic]]))
            rows.list[[topic]] <- list()
          if (is.null(rows.list_de[[topic]]))
            rows.list_de[[topic]] <- list()
          
          #print and enter
          temp <- mir_fam_go_list[[nam]]
          print(temp[temp$GO.ID %in% pertinent_ids,])
          message(nam)
          nums <- readline("ENTER ROWS: ")
          if (nchar(nums) > 0)
            rows.list[[topic]][[nam]] <-
            as.numeric(unlist(strsplit(nums, ",")))
          
          temp <- mir_fam_go_list_de[[nam]]
          print(temp[temp$GO.ID %in% pertinent_ids, ])
          message(nam)
          numsDE <- readline("ENTER ROWS (DE): ")
          if (nchar(nums) > 0)
            rows.list_de[[topic]][[nam]] <-
            as.numeric(unlist(strsplit(numsDE, ",")))
          count <- count + 1
        } else {
          count <- count + 1
        }
      }
    }
  }
}

save(topics, rows.list, rows.list_de, file = "working_data/go_curation.RData")
rows.list

ok <- names(mir_fam_go_list_de)
for(nam in ok){
  idx <- unlist(lapply(rows.list, function(x) return(x[nam])))
  idx_de <- unlist(lapply(rows.list_de, function(x) return(x[nam])))
  
  temp <- mir_fam_go_list[[nam]]
  if(length(idx)>0){
    print(temp[idx,])
    temp <- temp[-idx,]
  }
  print(temp[temp$GO.ID %in% pertinent_ids,])
  message(nam)
  topic <- readline("WHICH TOPIC")
  rows.list[[topic]][[nam]]
  
  temp <- mir_fam_go_list_de[[nam]]
  if(length(idx_de)>0){
    print(temp[idx_de,])
    temp <- temp[-idx_de,]
  }
  print(temp[temp$GO.ID %in% pertinent_ids,])
  topic <- readline("WHICH TOPIC")
  rows.list_de[[topic]][[nam]]
  
  readline("OK?")
  ok <- ok[-which(ok == nam)]
}


save(topics, rows.list, rows.list_de, file = "working_data/go_curation.RData")

#>heatmap####
mir_fam_go_list <- readRDS(file = "working_data/mir_family_GO_enrich_list.rds")
mir_fam_go_list_de <- readRDS(file = "working_data/mir_family_GO_enrich_list_de.rds")
df2 <- readRDS("working_data/de_enriched_families.rds")
load(file = "working_data/go_curation.RData")
topics_curation <- matrix(NA, nrow = nrow(df2), ncol = length(topics))
rownames(topics_curation) <- df2$name
colnames(topics_curation) <- sort(topics)
topics_curation <- topics_curation[, 1:(ncol(topics_curation)-1)]
topics_curation_de <- topics_curation

for(col in colnames(topics_curation)){
  for(row in rownames(topics_curation)){
    tab <- mir_fam_go_list[[row]]
    idx <- rows.list[[col]][[row]]
    idx_de <- rows.list_de[[col]][[row]]
    topics_curation[row, col] <- length(idx)
    topics_curation_de[row, col] <- length(idx_de)
  }
}

pheatmap::pheatmap(topics_curation)
pheatmap::pheatmap(topics_curation_de)

#both
topics_curation_b <- topics_curation
for(col in colnames(topics_curation_b)){
  for(row in rownames(topics_curation_b)){
    tab <- mir_fam_go_list[[row]]
    tab_de <- mir_fam_go_list_de[[row]]
    idx <- rows.list[[col]][[row]]
    idx_de <- rows.list_de[[col]][[row]]
    ids <- unique(c(tab$GO.ID[idx], tab_de$GO.ID[idx_de]))
    
    topics_curation_b[row, col] <- length(ids)
  }
}
pheatmap::pheatmap(t(topics_curation_b))

max(topics_curation_b)
col <- viridis::viridis(8)[c(1, 3:8)]
pheatmap::pheatmap(t(topics_curation_b),
                   cellwidth = 10,
                   cellheight = 20,
                   color = col,
                   filename = "img/go_curated_heatmap.pdf")


#>count####
sort(colSums(topics_curation_b))

#print####
for(fam in names(mir_fam_go_list)){
  message(fam)
  df <- mir_fam_go_list[[fam]]
  df_de <- mir_fam_go_list_de[[fam]]
  for(topic in names(rows.list)){
    message(topic)
    idx <- rows.list[[topic]][[fam]]
    if(length(idx)>0) {
      message("all")
      print(df[idx,])
    }
    idx_de <- rows.list_de[[topic]][[fam]]
    if(length(idx_de)>0) {
      message("de")
      print(df_de[idx_de,])
    }
  }
  readline("ENTER")
}