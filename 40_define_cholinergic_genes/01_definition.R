#DEFINE CHOLINERGIC SYSTEM FROM TINS REVIEW AND TARGETING TFs FROM MARBACH####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#libs####
library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")
library(biomaRt)
library(plyr)

# #biomart instance####
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(as.character(datasets$dataset[grep("Human", datasets$description)]), mart = ensembl)


#GENE SET DEFINITION####
#>>core genes####
#dataframe with name (SYMBOL), entrez, ensembl ids
core.genes.entrez <- c(1103,
                43,590,
                8292, #colq
                60482,6572,145270, #core
                
                1134,1135,1136,1137,1138,8973,1139,55584,57053,1140,1141,1142,1143,
                1144,1145,1146,1128,1129,1131,1132,1133, #receptor
                
                4804,#p75ntr = ngfr
                4914,#trka = ntrk1
                4915,#trkb #neurotrophic receptors
                1271, #cntfr
                3572, #gp130 = il6st
                3977,#lifr
                3570,#il6r #neurokine receptors
                4803,#ngf
                627, #bdnf 
                497258, #bdnf-as #neurotrophins
                1270, #cntf
                3976,#lif
                3569,#il-6 #neurokines
                2864, 2867, 2865, 338557,#ffars
                7080, #nkx2-1
                26468, #lhx6
                431707, #lhx8
                2636, #gbx1
                2637, #gbx1
                3670, #isl1 #development
                
                3187, 3188, #hnrnp h 1, 2
                1861, #torsin a
                1859, #dyrk1a
                1375, #cpt1b
                56261, #gpcpd1
                5130, #pcyt1a
                6558, #slc12a2
                47, #acly
                1119, #chka
                1120) #chkb #extended

core.genes <- data.frame()
for(i in 1:length(core.genes.entrez)){
  core.genes <- rbind(core.genes, cypher(graph, paste("MATCH (g:GENE) WHERE g.id =", core.genes.entrez[i],
                                         "RETURN g.id as entrez, g.name as gene_symbol, g.ensg as ensg")))
}

#manual corrections, 2 ensgs: A7 and B1
core.genes$ensg[core.genes$gene_symbol == "CHRNA7"] <- "ENSG00000175344"
core.genes$ensg[core.genes$gene_symbol == "CHRNB1"] <- "ENSG00000170175"
core.genes <- rbind(core.genes, c(89832, "CHRFAM7A", "ENSG00000166664"))

rownames(core.genes) <- seq(1, length(core.genes$entrez))
core.genes$group <- factor(c(rep("core", 7), 
                             rep("receptor", 21), 
                             rep("nt_receptor", 3), 
                             rep("nk_receptor", 4), 
                             rep("nt", 3), 
                             rep("nk", 3), 
                             rep("ffar", 4), 
                             rep("development", 6), 
                             rep("extended", 12)))

cholinergic_core_genes <- core.genes

#>>TFs####
#dataframe with name (SYMBOL), entrez, ensembl ids
#marbach from graph
if(!file.exists("./working_data/cholinergic_tfs_marbach_raw.rds")){
  tfs.marbach <-
    cypher(
      graph,
      "MATCH (tf:GENE {tf: TRUE})-[r]->(g:GENE)
      WHERE g.id in {ids}
      WITH tf, g, r, type(r) as type
      RETURN tf.name as tf_name, tf.id as tf_id, g.name as target_name, g.id as target_id, type as tissue, r.tfa as tfa",
      ids = as.numeric(cholinergic_core_genes$entrez)
    )
  saveRDS(tfs.marbach, file = "./working_data/cholinergic_tfs_marbach_raw.rds")
  } else {
    tfs.marbach <-
      readRDS("./working_data/cholinergic_tfs_marbach_raw.rds"
      )
  }
#filter marbach list by tfa ####
length(unique(tfs.marbach$tf_name))
hist(tfs.marbach$tfa)
quantile(tfs.marbach$tfa, seq(0,1,.01))
threshold <- .2
tfs.marbach <- tfs.marbach[tfs.marbach$tfa>threshold,]
length(unique(tfs.marbach$tf_name))
entrez2ensg <- getBM(attributes = c("entrezgene", "ensembl_gene_id"), filters = "entrezgene", 
                     values = c(tfs.marbach$tf_id, tfs.marbach$target_id), mart = ensembl)
tfs.marbach$tf_ensg <- entrez2ensg$ensembl_gene_id[match(tfs.marbach$tf_id, entrez2ensg$entrezgene)]
tfs.marbach$target_ensg <- entrez2ensg$ensembl_gene_id[match(tfs.marbach$target_id, entrez2ensg$entrezgene)]
head(tfs.marbach)
length(unique(tfs.marbach$tf_name))

#supplement tf list manually####
tfs.manual <- c("ISL1")
tfs.manual <- cypher(graph, "MATCH (g:GENE)
                     WHERE g.name in {names}
                     RETURN g.name as tf_name, g.id as tf_id", 
                     names = tfs.manual)
entrez2ensg <- getBM(attributes = c("entrezgene", "ensembl_gene_id"), filters = "entrezgene", values = tfs.manual$tf_id, mart = ensembl)
tfs.manual$tf_ensg <- entrez2ensg$ensembl_gene_id[match(tfs.manual$tf_id, entrez2ensg$entrezgene)]
tfs.manual

#supplement tf TRRUST####
trrust <- read.csv("raw_data/trrust_rawdata.human.tsv", sep = "\t", header = F)
head(trrust)
colnames(trrust) <- c("tf_name", "target_name", "effect", "PMID")
tfs.trrust <- data.frame(trrust[trrust$target_name %in% core.genes$gene_symbol,1:2])
trrust.ids <- cypher(graph, "MATCH (g:GENE)
                     WHERE g.name in {names}
                     RETURN g.name as gene_symbol, g.id as entrez", 
                     names = c(tfs.trrust$tf_name, tfs.trrust$target_name))
trrust.ids
tfs.trrust$tf_id <- trrust.ids$entrez[match(tfs.trrust$tf_name, trrust.ids$gene_symbol)]
tfs.trrust$target_id <- trrust.ids$entrez[match(tfs.trrust$target_name, trrust.ids$gene_symbol)]
tfs.trrust <- tfs.trrust[,c(1,3,2,4)]
tfs.trrust

entrez2ensg <- getBM(attributes = c("entrezgene", "ensembl_gene_id"), filters = "entrezgene", values = trrust.ids$entrez, mart = ensembl)
tfs.trrust$tf_ensg <- entrez2ensg$ensembl_gene_id[match(tfs.trrust$tf_id, entrez2ensg$entrezgene)]
tfs.trrust$target_ensg <- entrez2ensg$ensembl_gene_id[match(tfs.trrust$target_id, entrez2ensg$entrezgene)]
head(tfs.trrust)

#combine TF dataframes for targeting ####
head(tfs.marbach)
head(tfs.trrust)
head(tfs.manual)
tfs.comb <- rbind.fill(tfs.marbach, tfs.trrust, tfs.manual)

cholinergic_tfs <- tfs.comb
saveRDS(cholinergic_tfs, file = "working_data/cholinergic_tfs.rds")

#>>save TABLES####
save(cholinergic_core_genes, cholinergic_tfs, file = "./working_data/ecs_tables.RData")
