## SOURCE SINGLE DATASET DGE ####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)
options(stringsAsFactors = F)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

#set keys (ensg)####
cholinergic_genes <- readRDS("working_data/cholinergic_genes.rds")
cholinergic_tfs <- readRDS("working_data/cholinergic_tfs.rds")
tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")
mean(tfs$tf.name %in% cholinergic_tfs$tf_name)
unique(cholinergic_tfs$tf_ensg)
tfs_ensg <- cypher(graph, "MATCH (g:GENE) WHERE g.name IN {tfs} RETURN g.name, g.ensg", tfs = tfs$tf.name)

keys <- c(cholinergic_genes$ensg[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")], tfs_ensg$g.ensg)
saveRDS(keys, "working_data/subset_keys_ensg.rds")

length(unique(cholinergic_genes$ensg[!cholinergic_genes$group %in% c("extended", "ffar", "splicing")]))
nrow(tfs)

{
source("60_subset_analysis/01_chen/chen_read_affy.R")
source("60_subset_analysis/02_lanz/lanz_read_affy.R")
source("60_subset_analysis/03_iwamoto/iwamoto_read_affy.R")
source("60_subset_analysis/04_maycox/maycox_read_affy.R")
source("60_subset_analysis/05_narayan/narayan_read_affy.R")
source("60_subset_analysis/06_bahn/bahn_read_affy.R")
}
