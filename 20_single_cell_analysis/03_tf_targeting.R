#TF TARGETING SINGLE CELL SEQ DATA####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

genes <- cypher(graph, "MATCH (g:GENE {species: 'HSA'}) RETURN g.name")$g.name

darmanis <- readRDS("./raw_data/darmanis_chol_exprs.rds")
habib <- readRDS("./raw_data/habib_chol_exprs.rds")
tasic <- readRDS("./raw_data/tasic_chol_exprs.rds")
zeisel <- readRDS("./raw_data/zeisel_chol_exprs.rds")

t_query <- "MATCH (t:TISSUE {cell_class: 'nervous_cells'}) return t.name"
tissues <- toupper(cypher(graph, t_query)$t.name)
adult <- tissues[grep("ADULT", tissues)]
adult <- c(adult, "MANUAL_TF")

if(!file.exists("./working_data/single_cell/all_tfs_cns.rds")) {
  human_tf_frame <- data.frame()
  for (tissue in adult) {
    print(tissue)
    query <-
      paste0(
        "MATCH (tf:GENE {tf: TRUE})-[r:",
        tissue,
        "]-(g:GENE) ",
        "RETURN tf.name, type(r) as tissue, g.name, r.tfa"
      )
    human_tf_frame <- rbind(human_tf_frame,
                            cypher(graph, query))
  }
  head(human_tf_frame)
  saveRDS(human_tf_frame,
          "./working_data/single_cell/all_tfs_cns.rds")
} else {
  human_tf_frame <-
    readRDS("./working_data/single_cell/all_tfs_cns.rds")
}
#remove marker genes
marker <- rownames(darmanis)[36:42]
rowMeans(darmanis[1:35,1:9])
darmanis <- darmanis[1:35,]
rowMeans(habib[8:42,1:9])
habib <- habib[!habib$name %in% marker,]
rowMeans(tasic[!rownames(tasic) %in% marker,1:4])
tasic <- tasic[!tasic$name %in% marker,]
rowMeans(zeisel[!rownames(zeisel) %in% marker,1:4])
zeisel <- zeisel[!zeisel$name %in% marker,]

#more than one per row
darmanis <- darmanis[colSums(apply(darmanis[,1:9]>0, 1, as.numeric))>1,]
habib <- habib[colSums(apply(habib[,1:9]>0, 1, as.numeric))>1,]
tasic <- tasic[colSums(apply(tasic[,1:4]>0, 1, as.numeric))>1,]
zeisel <- zeisel[colSums(apply(zeisel[,1:4]>0, 1, as.numeric))>1,]

human <- rownames(darmanis)

mouse <- unique(c(habib$name, tasic$name, zeisel$name))

#human genes####
if(!file.exists("./working_data/single_cell/human_permutation_tf.rds")) {
  #tf targeting
  idx <- which(human_tf_frame$g.name %in% human)
  chol_tf_frame <- human_tf_frame[idx,]
  head(chol_tf_frame)
  nrow(chol_tf_frame)
  chol_tf_frame <- chol_tf_frame[chol_tf_frame$r.tfa>quantile(chol_tf_frame$r.tfa, seq(0,1,.01)[60]),] #min 50%
  chol_tf_frame <- chol_tf_frame[order(chol_tf_frame$r.tfa, decreasing = T),]
  chol_tf_frame <- chol_tf_frame[!duplicated(paste0(chol_tf_frame$tf.name, chol_tf_frame$g.name)),]
  human.count <- data.frame(dplyr::count(chol_tf_frame, tf.name, sort = T))
  head(human.count, 100)
  
  #permutation
  human.count$p.val <- NA
} else {
  human.count <-
    readRDS(file = "./working_data/single_cell/human_permutation_tf.rds")
}

if(ncol(human.count)<20000) {
  for (i in (ncol(human.count) + 1):20000) {
    print(i)
    perm_genes <- sample(genes, length(human))
    idx <- which(human_tf_frame$g.name %in% perm_genes)
    perm_frame <- human_tf_frame[idx,]
    perm_frame <-
      perm_frame[perm_frame$r.tfa > quantile(perm_frame$r.tfa, seq(0, 1, .1)[6]), ] #min 50%
    perm_frame <-
      perm_frame[order(perm_frame$r.tfa, decreasing = T), ]
    perm_frame <-
      perm_frame[!duplicated(paste0(perm_frame$tf.name, perm_frame$g.name)), ]
    targets.count <- data.frame(dplyr::count(perm_frame, tf.name))
    idx <- match(human.count$tf.name, targets.count$tf.name)
    human.count[, i] <- targets.count$n[idx]
    saveRDS(human.count, file = "./working_data/single_cell/human_permutation_tf.rds")
  }
} else {
  human.count <- readRDS(file = "./working_data/single_cell/human_permutation_tf.rds")
}


human.count[is.na(human.count)] <- 0
for (tf in human.count$tf.name) {
  null <-
    unlist(human.count[human.count$tf.name == tf, 4:ncol(human.count)])
  null <- sort(null)
  null <-
    cbind(null, data.frame(prob = 1 - abs(2 * seq(1:length(
      null
    )) / length(null) - 1)))
  count <- human.count$n[human.count$tf.name == tf]
  p <- null$prob[findInterval(count, null$null)]
  if (length(p) > 0)
    human.count$p.val[human.count$tf.name == tf] <- p
  else
    human.count$p.val[human.count$tf.name == tf] <- 0
}

human.count <- human.count[order(human.count$n, decreasing = T),]
rownames(human.count) <- seq(nrow(human.count))

tfs <- human.count[human.count$p.val<.05,c("tf.name", "n", "p.val")]
tfs
nrow(tfs)
tfs[tfs$p.val<.001,]

saveRDS(tfs, file = "./working_data/single_cell/significant_tfs.rds")

#SUPPLEMENTAL FIGURE: TF activity distributions####
t_query <- "MATCH (t:TISSUE {cell_class: 'nervous_cells'}) return t.name"
tissues <- toupper(cypher(graph, t_query)$t.name)
adult <- tissues[grep("ADULT", tissues)]
human_tf_frame <- data.frame()
for (tissue in adult) {
  print(tissue)
  query <-
    paste0(
      "MATCH (tf:GENE {tf: TRUE})-[r:",
      tissue,
      "]-(g:GENE) ",
      "RETURN tf.name, type(r) as tissue, g.name, r.tfa"
    )
  human_tf_frame <- rbind(human_tf_frame,
                          cypher(graph, query))
}
head(human_tf_frame)
saveRDS(human_tf_frame, "./working_data/single_cell/all_tfs_cns.rds")
sample_frame <- human_tf_frame#[sample(1:nrow(human_tf_frame), 1000000),]
ggplot(sample_frame[sample_frame$r.tfa > quantile(sample_frame$r.tfa, seq(0, 1, .01))[100],], 
       aes(r.tfa, color = tissue)) + geom_density() + facet_wrap("tissue") + 
  xlab(label = "Transcriptional activity") + ylab(label = "Density") +
  theme(legend.position = "none")
ggsave("./img/sup_fig_tfa_densities_top_1_perc.svg", width = 20, height = 15)
