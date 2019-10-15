#PREDICTED MIRS FROM CHOLINERGIC SINGLE CELLS####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(ggplot2)
library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")

long_term <- rownames(mir.matrix)[rowSums(mir.matrix[, 1:8]) != 0]

#single cell mirs####
primate_specific_mirs <- readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <- readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
tf_mirs <- readRDS("./working_data/single_cell/significant_tf_mirs.rds")

#primate
idx <- lapply(primate_specific_mirs$m.name, function(x) grep(paste0(x, "$"), long_term))
ps_found <- unique(unlist(lapply(idx, function(x) long_term[x])))
primate_specific_mirs[primate_specific_mirs$m.name %in% ps_found,]
ps_p <- nrow(primate_specific_mirs)
ps_f <- length(ps_found)
length(ps_found)/length(long_term)

ps <- data.frame(set = "primate", predicted = ps_p, de = ps_f)

#conserved
idx <- lapply(non_ps_mirs$m.name, function(x) grep(paste0(x, "$"), long_term))
nps_found <- unique(unlist(lapply(idx, function(x) long_term[x])))
non_ps_mirs[non_ps_mirs$m.name %in% nps_found & non_ps_mirs$p.val < .001,]
co_p <- nrow(non_ps_mirs)
co_f <- length(nps_found)
length(nps_found)/length(long_term)

co <- data.frame(set = "conserved", predicted = co_p, de = co_f)

#tf targeting
idx <- lapply(tf_mirs$m.name, function(x) grep(paste0(x, "$"), long_term))
tfm_found <- unique(unlist(lapply(idx, function(x) long_term[x])))
tf_mirs[tf_mirs$m.name %in% tfm_found & tf_mirs$p.val < .001,]
tf_p <- nrow(tf_mirs)
tf_f <- length(tfm_found)
length(tfm_found)/length(long_term)

tf <- data.frame(set = "TF", predicted = tf_p, de = tf_f)

ps_found[ps_found %in% tfm_found] #primate-specific TF mirs
nps_found[nps_found %in% tfm_found] #conserved TF mirs

#save
saveRDS(ps_found, file = "./working_data/mirna_seq/predicted_mirs_primate_overlap.rds")
saveRDS(nps_found, file = "./working_data/mirna_seq/predicted_mirs_conserved_overlap.rds")
saveRDS(tfm_found, file = "./working_data/mirna_seq/predicted_mirs_tf_mir_overlap.rds")

#bar graph
plot <- rbind(ps, co, tf)
plot$sum <- plot$predicted
plot$predicted <- (plot$predicted - plot$de) / plot$sum * 100
plot$de <- plot$de / plot$sum * 100
plot <- reshape2::melt(plot, id.vars = c("set", "sum"))
plot$set <- factor(plot$set, levels = c("conserved", "primate", "TF"))
ggplot(plot, aes(set, value, fill = variable)) + geom_bar(stat = "identity", position = "stack") +
  annotate("text", label = plot$sum[match(levels(plot$set), plot$set)], x = c(1,2,3), y = 80, size = 7, col = "black")
ggsave("img/prediction_de_overlap.pdf", width = 4, height = 5)

#mark in figure 3f####
cholinergic_genes <- readRDS("working_data/cholinergic_genes.rds")
miFam <- readRDS("working_data/mir_family_table_extended.rds")
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
ps_found

mirs <- unique(c(ps_found, nps_found, tfm_found))
genes <- cholinergic_genes$gene_symbol[!cholinergic_genes$group %in% c("extended", "ffar", "splicing", "BMP", "nt", "nt_receptor")]
targets <- get.targeting(mirs, genes)
targets <- unique(targets[, c("m.name", "g.name")])
nrow(targets) 
data.frame(dplyr::count(targets, m.name, sort = T))
targets[targets$m.name == "hsa-miR-1185-5p",]
targets[targets$m.name == "hsa-miR-320d",]
targets[targets$m.name == "hsa-miR-199a-5p",]
targets[targets$m.name == "hsa-miR-199b-5p",]
targets[targets$m.name == "hsa-miR-27b-5p",]
targets[targets$m.name == "hsa-miR-125b-5p",]
targets[targets$m.name == "hsa-let-7f-2-3p",]
 