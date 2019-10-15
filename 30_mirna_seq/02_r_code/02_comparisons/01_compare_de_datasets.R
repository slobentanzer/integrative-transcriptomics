#COMPARE DE DATASETS####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

options(stringsAsFactors = F)

library(RNeo4j)
library(venn)
library(reshape2)
library(RColorBrewer)
library(dplyr)
# library(pheatmap)
# library(rgexf)
# graph <- startGraph("http://localhost:7474/db/data/")
# cholinergic_genes <- readRDS("working_data/cholinergic_genes.rds")

#new sequencing####
# load DE
lnames <- load("working_data/mirna_seq/LA2_LA5_DESeq_results_new_seq.RData")
de_list <- vector(mode = "list", length = length(lnames))
names(de_list) <- lnames
for(nam in lnames){
  de_list[[nam]] <- get(nam)
}

all.mirs <- sort(unique(unlist(lapply(de_list, function(x) x$name))))

# matrix of timepoints
mir.matrix <- matrix(ncol = length(lnames), nrow = length(all.mirs), rep(0))
rownames(mir.matrix) <- all.mirs
mir.matrix <- data.frame(mir.matrix)
colnames(mir.matrix) <- lnames
#countchange
for(nam in colnames(mir.matrix)) {
  de <- de_list[[nam]]
  mir.matrix[, nam] <- de$countChange[match(rownames(mir.matrix), de$name[de$padj<.1])]
}
mir.matrix[is.na(mir.matrix)] <- 0

saveRDS(mir.matrix, file = "working_data/mir_de_matrix_countchange.rds")

#bar graph####
mm.melt <- melt(cbind(mir.matrix[,1:8], data.frame(name = rownames(mir.matrix))))
mm.melt <- mm.melt[mm.melt$value != 0,]
mm.melt %>% group_by(variable) %>% dplyr::count()
mm.melt$time <- unlist(lapply(strsplit(as.character(mm.melt$variable), "_"), "[", 2))
mm.melt$time <- factor(mm.melt$time, levels = c("min30", "min60", "day2", "day4"))
# mm.melt <- mm.melt[grep("day", mm.melt$time),]
mm.melt$cell <- unlist(lapply(strsplit(as.character(mm.melt$variable), "_"), "[", 1))
inter <- intersect(mm.melt$name[mm.melt$cell == "la2"], mm.melt$name[mm.melt$cell == "la5"])
mm.sep <- mm.melt[mm.melt$name %in% inter,]
mm.sep <- mm.sep[order(mm.sep$name),]
for(mir in unique(mm.sep$name)){
  temp <- mm.sep[mm.sep$name == mir,]
  # for(time in unique(temp$time)){
  df <- temp#[temp$time == time,]
  if(nrow(df) == 2){
    if(df$value[1] * df$value[2] < 0){
      df$cell <- "opposite"
    } else {
      df$cell <- "similar"
    }
  }
  if(nrow(df) > 2){
    if(any(df$value>0) & any(df$value<0)){
      df$cell <- "opposite"
    } else {
      df$cell <- "similar"
    }
  }
  temp <- df#[temp$time == time,] <- df
  # }
  mm.sep[mm.sep$name == mir,] <- temp
}

mm.sep %>% group_by_at(vars(one_of(c("time", "cell")))) %>% dplyr::count()

mm.sep <- mm.sep[order(mm.sep$value, decreasing = T),]
rem <- which(duplicated(mm.sep[, c(1,4)]))
mm.count <- rbind(mm.sep[-rem, ], mm.melt[!mm.melt$name %in% inter,])
mm.count <- mm.count[order(mm.count$name),]
which(duplicated(mm.count[, c(1,4)]))

mm.count$cell <- factor(mm.count$cell, levels = c("opposite", "similar", "la2", "la5"))
plot <- mm.count
plot$cell <- as.character(plot$cell)
plot$cell[plot$cell %in% c("opposite", "similar")] <- "both"
plot$cell <- factor(plot$cell, levels = c("la2", "both", "la5"))

ggplot(plot, aes(cell, fill = cell)) + geom_bar(position=position_dodge(width=0.9)) +
  facet_grid(.~time, scale="free_x", space="free")
ggsave("img/DE_mirs_per_timepoint_bar.svg", width = 5, height = 2)

#same, but for further reference of only long-term####
mm.melt <- melt(cbind(mir.matrix[,1:8], data.frame(name = rownames(mir.matrix))))
mm.melt <- mm.melt[mm.melt$value != 0,]
mm.melt %>% group_by(variable) %>% dplyr::count()
mm.melt$time <- unlist(lapply(strsplit(as.character(mm.melt$variable), "_"), "[", 2))
mm.melt$time <- factor(mm.melt$time, levels = c("min30", "min60", "day2", "day4"))
mm.melt <- mm.melt[grep("day", mm.melt$time),]
mm.melt$cell <- unlist(lapply(strsplit(as.character(mm.melt$variable), "_"), "[", 1))
inter <- intersect(mm.melt$name[mm.melt$cell == "la2"], mm.melt$name[mm.melt$cell == "la5"])
mm.sep <- mm.melt[mm.melt$name %in% inter,]
mm.sep <- mm.sep[order(mm.sep$name),]
for(mir in unique(mm.sep$name)){
  temp <- mm.sep[mm.sep$name == mir,]
  # for(time in unique(temp$time)){
  df <- temp#[temp$time == time,]
  if(nrow(df) == 2){
    if(df$value[1] * df$value[2] < 0){
      df$cell <- "opposite"
    } else {
      df$cell <- "similar"
    }
  }
  if(nrow(df) > 2){
    if(any(df$value>0) & any(df$value<0)){
      df$cell <- "opposite"
    } else {
      df$cell <- "similar"
    }
  }
  temp <- df#[temp$time == time,] <- df
  # }
  mm.sep[mm.sep$name == mir,] <- temp
}

mm.sep %>% group_by_at(vars(one_of(c("time", "cell")))) %>% dplyr::count()

mm.sep <- mm.sep[order(mm.sep$value, decreasing = T),]
rem <- which(duplicated(mm.sep[, c(1,4)]))
mm.count <- rbind(mm.sep[-rem, ], mm.melt[!mm.melt$name %in% inter,])
mm.count <- mm.count[order(mm.count$name),]
which(duplicated(mm.count[, c(1,4)]))

mm.count$cell <- factor(mm.count$cell, levels = c("opposite", "similar", "la2", "la5"))
plot <- mm.count
plot$cell <- as.character(plot$cell)
plot$cell[plot$cell %in% c("opposite", "similar")] <- "both"
plot$cell <- factor(plot$cell, levels = c("la2", "both", "la5"))

ggplot(plot, aes(cell, fill = cell)) + geom_bar(position=position_dodge(width=0.9)) +
  facet_grid(.~time, scale="free_x", space="free")

#count all unique DE mirs####
mm.count <- mm.count[order(abs(mm.count$value), decreasing = T),]
mm.count <- mm.count[order(mm.count$cell),]
mm.count <- mm.count[order(mm.count$name),]
mm.count.unique <- mm.count[!duplicated(mm.count$name),]

mm.count %>% group_by(cell) %>% dplyr::count()
mm.count.unique %>% group_by(cell) %>% dplyr::count()
nrow(mm.count.unique)

nrow(mir.matrix[rowSums(mir.matrix[,1:8])!=0,]) == nrow(mm.count.unique)

mm.count %>% group_by_at(vars(one_of(c("time", "cell")))) %>% dplyr::count()



#save list####
head(mm.sep)
mm.sep <- rbind(mm.sep, mm.melt[!mm.melt$name %in% inter,])
nrow(mm.sep) == nrow(mm.melt)

saveRDS(mm.sep, "working_data/de_counts_per_mir_condition_cell_line.rds")
saveRDS(mm.count, "working_data/de_counts_per_mir_unique.rds")

#countchange comparison####
la2 <- mir.matrix[, grep("la2_", colnames(mir.matrix))]
la5 <- mir.matrix[, grep("^la5_", colnames(mir.matrix))]
all(rownames(la5) == rownames(la2))
plot(asinh(c(la5[, 1], la5[, 2])) ~ asinh(c(la2[, 1], la2[, 2]))) 

#correlation of overlap####
la2 <- c(la2[, 1], la2[, 2])
la5 <- c(la5[, 1], la5[, 2])

#same direction
idx_same <- which(la2 * la5 > 0)
plot(asinh(la2[idx_same]), asinh(la5[idx_same]))
cor.test(la2[idx_same], la5[idx_same], method = "spear")
#dif direction
idx_dif <- which(la2 * la5 < 0)
plot(asinh(la2[idx_dif]), asinh(la5[idx_dif]))
cor.test(la2[idx_dif], la5[idx_dif], method = "spear")
#single changed
idx_single <- which(la2 * la5 == 0)
plot(asinh(la2[idx_single]), asinh(la5[idx_single]))

svg("img/countChange_correlation.svg", height = 5, width = 5)
plot(asinh(la2[idx_same]), asinh(la5[idx_same]), pch = 24, xlim = c(-10, 15), ylim = c(-15, 20), col = "green")
points(asinh(la2[idx_single]), asinh(la5[idx_single]), pch = 22)
points(asinh(la2[idx_dif]), asinh(la5[idx_dif]), pch = 25, col = "yellow")
dev.off()


la2 <- la2[la2 != 0]
hist(la2)

la5 <- la5[la5 != 0]
hist(la5)

t.test(abs(la5), abs(la2), alternative = "g")

#only intersecting####
la2 <- mir.matrix[, grep("la2_", colnames(mir.matrix))]
la5 <- mir.matrix[, grep("^la5_", colnames(mir.matrix))]
la2 <- la2[rowSums(la2) != 0,]
la5 <- la5[rowSums(la5) != 0,]
la2 <- la2[rownames(la2) %in% rownames(la5),]
la5 <- la5[rownames(la5) %in% rownames(la2),]
all(rownames(la5) == rownames(la2))

plot(asinh(c(la5[, 1], la5[, 2])) ~ asinh(c(la2[, 1], la2[, 2]))) 
plot(asinh(c(la5[, 3], la5[, 4])) ~ asinh(c(la2[, 3], la2[, 4]))) 
plot(asinh(c(la5[, 1], la5[, 2], la5[, 3], la5[, 4])) ~ asinh(c(la2[, 1], la2[, 2], la2[, 3], la2[, 4]))) 

la2 <- la2[la2 != 0]
hist(la2)

la5 <- la5[la5 != 0]
hist(la5)

#binary
for(i in 1:length(lnames)) {
  de <- de_list[[i]]
  mir.matrix[which(rownames(mir.matrix) %in% de$name[de$padj < .1]), i] <- 1
}
mir.matrix

svg(paste0("img/", "la2_expression_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(3,4,1,2)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la2_4_timepoint_overlap <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(3,4,1,2)])==4, ]))

svg(paste0("img/", "la5_expression_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(7,8,5,6)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la5_4_timepoint_overlap <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(7,8,5,6)])==4, ]))


svg(paste0("img/", "la2_la5_long_term_overlap_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(1,2,5,6)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la2_la5_long_term_overlap <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(1,2,5,6)])==4, ]))

svg(paste0("img/", "la2_vs_la5_expression_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(11,12,9,10)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la2_vs_la5_4_timepoint_overlap <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(11,12,9,10)])==4, ]))

svg(paste0("img/", "la2_interaction_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(1,2,9,10)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la2_vs_la5_overlap_with_la2 <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(1,2,9,10)])==4, ]))

svg(paste0("img/", "la5_interaction_venn.svg"), width = 8, height = 8)
venn(data.frame(mir.matrix[c(5,6,9,10)]),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()
la2_vs_la5_overlap_with_la5 <- sort(rownames(mir.matrix[rowSums(mir.matrix[, c(5,6,9,10)])==4, ]))

intersect(la2_vs_la5_overlap_with_la2, la2_vs_la5_overlap_with_la5)

#validate la2 with old sequencing####
# load DE
lnames <- load("working_data/mirna_seq/LA2_LA5_DESeq_results_old_seq_validation.RData")
de_list <- vector(mode = "list", length = length(lnames))
names(de_list) <- lnames
for(nam in lnames){
  de_list[[nam]] <- get(nam)
}

# matrix of timepoints
mir.matrix_val <- matrix(ncol = length(lnames), nrow = length(all.mirs), rep(0))
rownames(mir.matrix_val) <- all.mirs
mir.matrix_val <- data.frame(mir.matrix_val)
colnames(mir.matrix_val) <- lnames

#binary
for(i in 1:length(lnames)) {
  de <- de_list[[i]]
  mir.matrix_val[which(rownames(mir.matrix_val) %in% de$name[de$padj < .1]), i] <- 1
}
mir.matrix_val
colnames(mir.matrix_val) <- paste0(colnames(mir.matrix_val), "_val")

svg(paste0("img/", "la2_validation_venn.svg"), width = 8, height = 8)
val_frame <- cbind(mir.matrix[, 1:2], mir.matrix_val[, 1:2])
venn(data.frame(val_frame),
     zcolor = brewer.pal(4, "Set1"), cexsn = .7, cexil = 1)
dev.off()

val_frame <- cbind(mir.matrix[, 1:4], mir.matrix_val[, 1:4])

val_perc <- length(which(apply(val_frame[, 5:8], 1, function(row) any(row > 0)))) /
  length(which(apply(val_frame[, 1:4], 1, function(row) any(row > 0)))) #77.8 % validated

valmirs <- names(which(apply(val_frame[, 5:8], 1, function(row) any(row > 0)))) #validated mirs
saveRDS(valmirs, "working_data/validated_de_mirs_la2.rds")

#difference to added old with batches####
# load DE
lnames <- load("working_data/mirna_seq/LA2_LA5_DESeq_results_all_seq_batches.RData")
de_list <- vector(mode = "list", length = length(lnames))
names(de_list) <- lnames
for(nam in lnames){
  de_list[[nam]] <- get(nam)
}

# matrix of timepoints
mir.matrix_new <- matrix(ncol = length(lnames), nrow = length(all.mirs), rep(0))
rownames(mir.matrix_new) <- all.mirs
mir.matrix_new <- data.frame(mir.matrix_new)
colnames(mir.matrix_new) <- lnames
# #countchange
# for(nam in colnames(mir.matrix_new)) {
#   de <- de_list[[nam]]
#   mir.matrix_new[, nam] <- de$countChange[match(rownames(mir.matrix_new), de$name[de$padj<.1])]
# }
# mir.matrix_new[is.na(mir.matrix_new)] <- 0

#binary
for(i in 1:length(lnames)) {
  de <- de_list[[i]]
  mir.matrix_new[which(rownames(mir.matrix_new) %in% de$name[de$padj < .1]), i] <- 1
}
mir.matrix_new

for(i in 1:ncol(mir.matrix)){
  temp <- data.frame(mir.matrix[i], mir.matrix_new[i])
  colnames(temp) <- paste(colnames(temp), c("NEW", "BATCHES"))
  venn(temp,
        zcolor = brewer.pal(3, "Set1"), cexsn = .7, cexil = 1)
  readline("ENTER")
}

la2.2_old <- rownames(mir.matrix)[mir.matrix$la2_day2 == 1]
la2.2_new <- rownames(mir.matrix_new)[mir.matrix_new$la2_day2 == 1]
la2.2_old[!la2.2_old %in% la2.2_new]
la2.2_new[!la2.2_new %in% la2.2_old]

la2.4_old <- rownames(mir.matrix)[mir.matrix$la2_day4 == 1]
la2.4_new <- rownames(mir.matrix_new)[mir.matrix_new$la2_day4 == 1]
la2.4_old[!la2.4_old %in% la2.4_new]
la2.4_new[!la2.4_new %in% la2.4_old]

