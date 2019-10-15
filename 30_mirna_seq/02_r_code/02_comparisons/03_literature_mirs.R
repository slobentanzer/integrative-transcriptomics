#LITERATURE MIRS IN SCZ AND BD####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(ggplot2)
library(venn)
library(RColorBrewer)
library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")


mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")
mir.matrix <- mir.matrix[rowSums(mir.matrix[,1:8]) != 0,]
nrow(mir.matrix) #DE mirs
nrow(mir.matrix[rowSums(mir.matrix[,1:4]) != 0,]) #DE mirs
nrow(mir.matrix[rowSums(mir.matrix[,5:8]) != 0,]) #DE mirs
long_term <- rownames(mir.matrix)

long_term <- unlist(lapply(strsplit(long_term, "-"), function(x) paste(x[1:3], collapse = "-")))
long_term <- gsub("r", "R", long_term)
long_term_lookup <- data.frame(shortname = long_term, name = rownames(mir.matrix)) 
long_term <- unique(long_term)

#literature mirs####
#SCZ####
{
bc2011mirs <- read.table("./raw_data/bc2011scz_mirs.txt", encoding = "UTF-8")
bc2011mirs <- bc2011mirs[,1]
bc2011mirs <- bc2011mirs[c(grep("miR‐", bc2011mirs, fixed = T), grep("let‐", bc2011mirs, fixed = T))]
bc2011mirs <- bc2011mirs[-grep("SNP", bc2011mirs, fixed = T)]
bc2011mirs <- gsub("(", "", bc2011mirs, fixed = T)
bc2011mirs <- gsub(")", "", bc2011mirs, fixed = T)
bc2011mirs <- gsub("*", "", bc2011mirs, fixed = T)
bc2011mirs <- gsub("‐", "-", bc2011mirs, fixed = T)
bc2011mirs <- gsub("-3p", "", bc2011mirs, fixed = T)
bc2011mirs <- gsub("-5p", "", bc2011mirs, fixed = T)
bc2011mirs[2] <- "miR-128"
bc2011mirs <- unique(bc2011mirs)
}
idx <- lapply(bc2011mirs, function(x) grep(paste0(x, "$"), long_term))
scz_found <- unique(unlist(lapply(idx, function(x) long_term[x])))
length(bc2011mirs)
length(long_term)
length(scz_found)
length(scz_found)/length(long_term)

#BD####
{
fries2018mirs <- scan("./raw_data/fries2018bd_mirs.txt", what = "")
fries2018mirs <- fries2018mirs[c(grep("miR-", fries2018mirs, fixed = T), grep("let-", fries2018mirs, fixed = T))]
fries2018mirs <- gsub("*", "", fries2018mirs, fixed = T)
fries2018mirs <- gsub(",", "", fries2018mirs, fixed = T)
fries2018mirs <- gsub(".", "", fries2018mirs, fixed = T)
fries2018mirs <- gsub(";", "", fries2018mirs, fixed = T)
fries2018mirs <- gsub("-3p", "", fries2018mirs, fixed = T)
fries2018mirs <- gsub("-5p", "", fries2018mirs, fixed = T)
fries2018mirs <- unique(fries2018mirs)
}

idx <- lapply(fries2018mirs, function(x) grep(paste0(x, "$"), long_term))
bd_found <- unique(unlist(lapply(idx, function(x) long_term[x])))
length(fries2018mirs)
length(bd_found)
length(bd_found)/length(long_term)

intersect_lit <- unique(fries2018mirs[fries2018mirs %in% bc2011mirs])
intersect_found <- bd_found[bd_found %in% scz_found]

#save####
save(bc2011mirs, fries2018mirs, long_term, file = "./working_data/mirna_seq/mir_disease_literature.RData")

#table####
both_mirs <- unique(c(bc2011mirs, fries2018mirs))
both_mirs <- both_mirs[order(both_mirs)]
both_mirs <- data.frame(name = both_mirs)

both_mirs$BD <- both_mirs$name %in% fries2018mirs
both_mirs$SCZ <- both_mirs$name %in% bc2011mirs

both_mirs$citation <- NA
for(i in 1:nrow(both_mirs)){
  BD <- as.numeric(both_mirs$BD[i])
  SCZ <- as.numeric(both_mirs$SCZ[i])
  if(BD+SCZ == 2){
    both_mirs$citation[i] <- "Fries et al., 2018 & Beveridge and Cairns, 2012"
  } else if (BD == 1) {
    both_mirs$citation[i] <- "Fries et al., 2018"
  } else {
    both_mirs$citation[i] <- "Beveridge and Cairns, 2012"
  }
}

both_mirs$name <- paste0("hsa-", both_mirs$name)

write.table(both_mirs, file = "./out/DataS4.csv", row.names = F, quote = T, sep = ";")

#predicted overlap####
primate_pred <- readRDS(file = "./working_data/mirna_seq/predicted_mirs_primate_overlap.rds")
conserved_pred <- readRDS(file = "./working_data/mirna_seq/predicted_mirs_conserved_overlap.rds")
tfm_pred <- readRDS(file = "./working_data/mirna_seq/predicted_mirs_tf_mir_overlap.rds")

primate_pred <- unlist(lapply(strsplit(primate_pred, "-"), function(x) paste(x[2:3], collapse = "-")))
conserved_pred <- unlist(lapply(strsplit(conserved_pred, "-"), function(x) paste(x[2:3], collapse = "-")))
tfm_pred <- unlist(lapply(strsplit(tfm_pred, "-"), function(x) paste(x[2:3], collapse = "-")))

primate_pred[primate_pred %in% bc2011mirs]
primate_pred[primate_pred %in% fries2018mirs]

conserved_pred[conserved_pred %in% bc2011mirs]
conserved_pred[conserved_pred %in% fries2018mirs]

tfm_pred[tfm_pred %in% bc2011mirs]
tfm_pred[tfm_pred %in% fries2018mirs]

#venn####
length(bc2011mirs[bc2011mirs %in% gsub("hsa-", "", long_term)])
length(fries2018mirs[fries2018mirs %in% gsub("hsa-", "", long_term)])
bc_fries_mirs <- fries2018mirs[fries2018mirs %in% bc2011mirs]
length(bc_fries_mirs)
length(bc_fries_mirs[bc_fries_mirs %in% gsub("hsa-", "", long_term)])
length(long_term)
length(bc2011mirs)
length(fries2018mirs)


allmirs <- unique(c(gsub("hsa-", "", long_term), bc2011mirs, fries2018mirs))

venn_frame <- matrix(ncol = 3, nrow = length(allmirs), rep(0))
rownames(venn_frame) <- allmirs
venn_frame <- data.frame(venn_frame)
colnames(venn_frame) <- c("DE", "SCZ", "BD")
#countchange
for(i in 1:3) {
  temp <- switch(i, gsub("hsa-", "", long_term), bc2011mirs, fries2018mirs)
  venn_frame[, i] <- as.numeric(rownames(venn_frame) %in% temp)
}
svg("img/literature_venn.svg")
venn(venn_frame,
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
venn_frame[rowSums(venn_frame) == 3,]
venn_frame[rowSums(venn_frame) == 2,]

saveRDS(venn_frame, file = "working_data/literature_de_overlap_frame.rds")
