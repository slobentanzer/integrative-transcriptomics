#MIR PROPERTIES, ENRICHMENT, DE####
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

#DIFFERENTIAL EXPRESSION####
miFam$de <- miFam$m.name %in% rownames(mir.matrix)
miFam$timepoints <- mir.matrix$de_timepoints[match(miFam$m.name, rownames(mir.matrix))]
miFam$timepoints[is.na(miFam$timepoints)] <- 0

#FAMILY HYPERGEOMETRIC####
sig_fams <- readRDS(file = "working_data/cholinergic_permutation_family_enrichment.rds")
chol_enrich <- names(which(unlist(lapply(sig_fams, function(x) x[[2]] < .05))))
miFam$neurokine_enriched <- miFam$family_name %in% chol_enrich

#SINGLE CELL ENRICHMENT####
primate_specific_mirs <-
  readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <-
  readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
tf_mirs <-
  readRDS(file = "./working_data/single_cell/significant_tf_mirs.rds")

miFam$single_cell_enriched <- miFam$m.name %in% c(primate_specific_mirs$m.name,
                                                  non_ps_mirs$m.name,
                                                  tf_mirs$m.name)
miFam$single_cell_tf_mir <- miFam$m.name %in% tf_mirs$m.name

#DE ENRICHMENT####
df2 <- readRDS("working_data/de_enriched_families.rds")
# miFam$de_enriched_both <- miFam$family_name %in% df2$name[df2$log_pval_la2>0 & df2$log_pval_la5>0]
miFam$de_enriched_female <- miFam$family_name %in% df2$name[df2$log_pval_la2>0]
miFam$de_enriched_male <- miFam$family_name %in% df2$name[df2$log_pval_la5>0]

#save####
saveRDS(miFam, file = "working_data/mir_family_table_extended.rds")

#VISUALIZE####
vis <- apply(miFam[, grep("neurokine", colnames(miFam)):ncol(miFam)], 2, as.numeric)
rownames(vis) <- miFam$m.name
vis <- vis[!duplicated(rownames(vis)),]
pheatmap::pheatmap(vis)
pheatmap::pheatmap(vis[rowSums(vis)>2,])

rownames(vis)[rowSums(vis)>2]
miFam.f <- miFam[miFam$m.name %in% rownames(vis)[rowSums(vis)>2],]
length(unique(miFam.f$m.name))
length(unique(miFam.f$family_name))
table(miFam.f$family_name)

#save supplemental Data S8####
S8 <- miFam.f[, c("m.name", "family_name", "family_accession", "de", 
                  "single_cell_enriched", "single_cell_tf_mir", "de_enriched_female", "de_enriched_male")]
write.table(S8, file = "out/DataS8.csv", row.names = F, quote = F, sep = ";")
