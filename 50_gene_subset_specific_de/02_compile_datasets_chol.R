#COMPILE SUBSET DIFFERENTIAL EXPRESSION FROM STUDIES####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

options(stringsAsFactors = F)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")


# #libs####
library(ggplot2)
cholinergic_genes <- readRDS("working_data/cholinergic_genes.rds")
cholinergic_genes <- cholinergic_genes[!cholinergic_genes$group %in% c("extended", "ffar", "splicing"),]
tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")


array.folder <- "./working_data/subset/dge/"
array.files <- paste0(array.folder, list.files(array.folder))
sex.files <- array.files[grep("sex", array.files)]
chol.sex <- sex.files[grep("chol", sex.files)]
mtar.sex <- sex.files[grep("mtar", sex.files)]

#>>CHOLINERGIC#####

#>SCZ####
scz.files <- array.files[grep("scz", array.files)]
scz.files <- scz.files[grep("chol", scz.files)]
scz.files <- scz.files[-grep("sex", scz.files)]

scz.names <- c("Chen", "Iwamoto", "Lanz", "Maycox", "Narayan")
scz_comp <- data.frame()
for (i in 1:length(scz.files)) {
  dgelist <- readRDS(scz.files[i])
  dgelist <- dgelist[grep("SCZ", names(dgelist))]
  dgelist <- do.call("rbind", dgelist)
  if (!is.null(dgelist)) {
    dgelist$set <- scz.names[i]
    dgelist <-
      dgelist[, c("name", "logFC", "adj.P.Val", "sex", "set")]
    
    scz_comp <- rbind(scz_comp, dgelist)
  }
}

ggplot(scz_comp[scz_comp$name %in% cholinergic_genes$gene_symbol,], aes(name, logFC, fill = set)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(
    angle = 60,
    hjust = 1,
    size = 12
  ))

#comp####
scz_comp <- scz_comp[order(scz_comp$adj.P.Val),]
scz_comp <- scz_comp[!is.na(scz_comp$name),]
scz_comp.agg <- scz_comp[!duplicated(scz_comp$name),]
table(scz_comp.agg$set)

ggplot(scz_comp.agg[scz_comp.agg$name %in% cholinergic_genes$gene_symbol,], aes(name, logFC, fill = set)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(
    angle = 60,
    hjust = 1,
    size = 12
  ))

#>BD####
bd.files <- array.files[grep("bd", array.files)]
bd.files <- bd.files[grep("chol", bd.files)]
bd.files <- bd.files[-grep("sex", bd.files)]

bd.names <- c("Bahn", "Iwamoto", "Lanz")
bd_comp <- data.frame()
for (i in 1:length(bd.files)) {
  dgelist <- readRDS(bd.files[i])
  dgelist <- dgelist[grep("BD", names(dgelist))]
  dgelist <- do.call("rbind", dgelist)
  if (!is.null(dgelist)) {
    colnames(dgelist)[colnames(dgelist) == "hgnc_symbol"] <- "name"
    dgelist$set <- bd.names[i]
    dgelist <-
      dgelist[, c("name", "logFC", "adj.P.Val", "sex", "set")]
    
    bd_comp <- rbind(bd_comp, dgelist)
  }
  
}

ggplot(bd_comp, aes(name, logFC, fill = set)) + geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

bd_comp <- bd_comp[order(bd_comp$adj.P.Val),]
bd_comp.agg <- bd_comp[!duplicated(bd_comp$name),]

ggplot(bd_comp.agg, aes(name, logFC, fill = set)) + geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

#>SEX####
#BD####
bd.sex <- chol.sex[grep("bd", chol.sex)]

bd.names <- c("Bahn", "Iwamoto", "Lanz")
bd_comp <- data.frame()
for (i in 1:length(bd.sex)) {
  dgelist <- readRDS(bd.sex[i])
  dgelist <- dgelist[grep("BD", names(dgelist))]
  dgelist <- do.call("rbind", dgelist)
  if (!is.null(dgelist)) {
    colnames(dgelist)[colnames(dgelist) == "hgnc_symbol"] <- "name"
    dgelist$set <- bd.names[i]
    dgelist <-
      dgelist[, c("name", "logFC", "adj.P.Val", "sex", "set")]
    
    bd_comp <- rbind(bd_comp, dgelist)
  }
}

bd_comp
bd_all <- rbind(bd_comp, bd_comp.agg)
length(unique(bd_all$name))

ggplot(bd_all, aes(sex, name)) + geom_count()

#SCZ####
scz.sex <- chol.sex[grep("scz", chol.sex)]

scz.names <- c("Chen", "Iwamoto", "Lanz", "Maycox", "Narayan")
scz_comp <- data.frame()
for (i in 1:length(scz.sex)) {
  dgelist <- readRDS(scz.sex[i])
  dgelist <- dgelist[grep("SCZ", names(dgelist))]
  dgelist <- do.call("rbind", dgelist)
  if (!is.null(dgelist)) {
    dgelist$set <- scz.names[i]
    dgelist <-
      dgelist[, c("name", "logFC", "adj.P.Val", "sex", "set")]
    
    scz_comp <- rbind(scz_comp, dgelist)
  }
}
scz_comp

scz_all <- rbind(scz_comp, scz_comp.agg)
length(unique(scz_all$name))

ggplot(scz_all, aes(sex, name)) + geom_count()

#compilation####
bd_all$condition <- "BD"
scz_all$condition <- "SCZ"

chol_all <- rbind(bd_all, scz_all)
chol_all$genes <- "cholinergic"

ggplot(chol_all, aes(sex, name)) + geom_count() + facet_wrap("condition")
length(unique(chol_all$name))


#SAVE CHOLINERGIC GENES####
saveRDS(chol_all, file = "./working_data/subset/chol_bd_scz.rds")
chol_all <- readRDS(file = "./working_data/subset/chol_bd_scz.rds")
  
# find TFs
gene_tf <-
  cypher(
    graph,
    "MATCH (g:GENE {species: 'HSA'}) WHERE g.name in {names} RETURN g.name, g.tf",
    names = chol_all$name
  )
gene_tf$g.tf[is.na(gene_tf$g.tf)] <- F

idx <- match(chol_all$name, gene_tf$g.name)
chol_all$biotype <- gene_tf$g.tf[idx]
chol_all[which(is.na(chol_all$biotype)), ]
chol_all$biotype[which(is.na(chol_all$biotype))] <- F

chol_all$biotype <- factor(chol_all$biotype)
levels(chol_all$biotype) <- c("gene", "TF")

chol_non <- data.frame(name = cholinergic_genes$gene_symbol[!cholinergic_genes$gene_symbol %in% chol_all$name], 
                       logFC = NA, 
                       adj.P.Val = NA,
                       sex = NA,
                       set = NA,
                       condition = NA,
                       biotype = "gene")
chol_tf <-
  cypher(
    graph,
    "MATCH (g:GENE {species: 'HSA'}) WHERE g.name in {names} RETURN g.name, g.tf",
    names = chol_non$name
  )
chol_non$biotype[chol_non$name %in% chol_tf$g.name[chol_tf$g.tf]] <- "TF"

tf_non <- data.frame(name = tfs$tf.name[!tfs$tf.name %in% chol_all$name], 
                     logFC = NA, 
                     adj.P.Val = NA,
                     sex = NA,
                     set = NA,
                     condition = NA,
                     biotype = "TF")

chol_all <- rbind(chol_all[, c("name","logFC","adj.P.Val","sex","set","condition","biotype")],
                  chol_non[, c("name","logFC","adj.P.Val","sex","set","condition","biotype")],
                  tf_non[, c("name","logFC","adj.P.Val","sex","set","condition","biotype")])

nrow(chol_all)
length(unique(chol_all$name))
chol_all <- chol_all[chol_all$name %in% c(cholinergic_genes$gene_symbol, tfs$tf.name),]

c(cholinergic_genes$gene_symbol, tfs$tf.name)[!c(cholinergic_genes$gene_symbol, tfs$tf.name) %in% chol_all$name]

#save supplementary Data S7####
write.table(chol_all, file = "out/DataS7.csv", row.names = F, quote = F, sep = ";")
