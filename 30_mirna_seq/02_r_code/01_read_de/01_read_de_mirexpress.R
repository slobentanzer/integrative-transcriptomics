#READ MIREXPRESS ALIGNMENT TO miRBase v21####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

options(stringsAsFactors = F)

library(DESeq2)


get.count.change <- function(BM, LF) {
  return((BM * 2^LF) - BM)
}

ancestors <- read.table(file = "./working_data/mirna_seq/mir_ancestors.csv", header = T, sep = "\t")

#READ####
#EXPERIMENT LA2_LA5####
mirexp.folder <- "./working_data/mirna_seq/mirexpress_alignment/LA2_LA5/results"
mirexp.samples <- sort(list.dirs(mirexp.folder, full.names = F, recursive = F))
mirexp.dirs <- sort(list.dirs(mirexp.folder, recursive = F))

n <- length(mirexp.dirs)
mirexp <- vector(mode = "list", length = n)
names(mirexp) <- mirexp.samples
for(i in 1:n){
  temp <- read.table(paste0(mirexp.dirs[i], "/miRNA_expression.txt"), sep = "\t")
  colnames(temp)[2] <- "count"
  temp$anc <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 1))
  temp$mat <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 2))
  temp$V1 <- NULL
  head(temp)
  temp <- temp[order(temp$count, decreasing = T),]
  temp <- temp[which(!duplicated(temp$mat)),]
  mirexp[[i]] <- temp
}

all_mirs <- sort(unique(unlist(lapply(mirexp, function(x) x$mat))))
counts <- data.frame(name = all_mirs)
for(i in 1:n) {
  idx <- match(counts$name, mirexp[[i]]$mat)
  counts[, (i+1)] <- mirexp[[i]]$count[idx]
}

rownames(counts) <- counts$name
counts$name <- NULL
counts[is.na(counts)] <- 0
names(counts) <- names(mirexp)
head(counts)

write.csv(counts, "out/mirexpress2_counts_raw.csv", row.names = T, quote = F)

## sample table ####
sample.table <- data.frame(sample = mirexp.samples)
rownames(sample.table) <- sample.table$sample
sample.table$sample.name <- unlist(lapply(strsplit((sample.table$sample), "-"), "[", 2))
sample.table$condition <- as.factor(grepl("CN", sample.table$sample.name))
levels(sample.table$condition) <- c("CON", "CNTF")

sample.table$time <- sample.table$sample.name
sample.table$time <- gsub("LA2CO", "", sample.table$time); sample.table$time <- gsub("LA5CO", "", sample.table$time); 
sample.table$time <- gsub("LA2CN", "", sample.table$time); sample.table$time <- gsub("LA5CN", "", sample.table$time); 
sample.table$time <- as.factor(sample.table$time)
levels(sample.table$time) <- c("DAY2", "MIN30", "DAY4", "MIN60")

sample.table$type <- as.factor(substr(sample.table$sample.name, 1, 3))

sample.table$con_time_type <- as.factor(apply(sample.table[, c("condition", "time", "type")], 1, paste, collapse = "_"))
sample.table$con_time <- as.factor(apply(sample.table[, c("condition", "time")], 1, paste, collapse = "_"))

dim(counts)
dim(sample.table)
all(colnames(counts) == rownames(sample.table))

#split sample tables per cell type and day for simple DE####
sample.table_la2 <- sample.table[sample.table$type == "LA2",]
sample.table_la2.2 <- sample.table_la2[sample.table_la2$time != "DAY4",]
sample.table_la2.2$con_time <- as.factor(as.character(sample.table_la2.2$con_time))
sample.table_la2.2$con_time <- relevel(sample.table_la2.2$con_time, "CON_DAY2")
sample.table_la2.4 <- sample.table_la2[sample.table_la2$time == "DAY4",]
sample.table_la2.4$con_time <- as.factor(as.character(sample.table_la2.4$con_time))
sample.table_la2.4$con_time <- relevel(sample.table_la2.4$con_time, "CON_DAY4")

sample.table_la5 <- sample.table[sample.table$type == "LA5",]
sample.table_la5.2 <- sample.table_la5[sample.table_la5$time != "DAY4",]
sample.table_la5.2$con_time <- as.factor(as.character(sample.table_la5.2$con_time))
sample.table_la5.2$con_time <- relevel(sample.table_la5.2$con_time, "CON_DAY2")
sample.table_la5.4 <- sample.table_la5[sample.table_la5$time == "DAY4",]
sample.table_la5.4$con_time <- as.factor(as.character(sample.table_la5.4$con_time))
sample.table_la5.4$con_time <- relevel(sample.table_la5.4$con_time, "CON_DAY4")


## DESeq SINGLE SETS DE ####
#>>LA2, DAY2####
se_la2.2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.2$sample]), 
                                 colData = sample.table_la2.2)
ddsSE_2.2 <- DESeqDataSet(se_la2.2, design = ~ con_time)
ddsSE_2.2

#>DESeq run ####
dds_2.2 <- estimateSizeFactors(ddsSE_2.2)
nc <- counts(dds_2.2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.2 <- dds_2.2[filter,]
dds_2.2 <- estimateDispersions(dds_2.2)
dds_2.2 <- nbinomWaldTest(dds_2.2, maxit=500)
resultsNames(dds_2.2)

#30 minutes vs control
la2_min30 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN30", "CON_DAY2"))
plotMA(la2_min30)
summary(la2_min30)

la2_min30$name <- rownames(la2_min30)
la2_min30$countChange <- get.count.change(la2_min30$baseMean, la2_min30$log2FoldChange)
la2_min30 <- la2_min30[order(la2_min30$padj, decreasing = F),]
la2_min30 <- la2_min30[!is.na(la2_min30$padj),]
data.frame(la2_min30)


#60 minutes vs control
la2_min60 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN60", "CON_DAY2"))
plotMA(la2_min60)
summary(la2_min60)

la2_min60$name <- rownames(la2_min60)
la2_min60$countChange <- get.count.change(la2_min60$baseMean, la2_min60$log2FoldChange)
la2_min60 <- la2_min60[order(la2_min60$padj, decreasing = F),]
la2_min60 <- la2_min60[!is.na(la2_min60$padj),]
data.frame(la2_min60)

#2 days vs control
la2_day2 <- results(dds_2.2, contrast = c("con_time", "CNTF_DAY2", "CON_DAY2"))
plotMA(la2_day2)
summary(la2_day2)

la2_day2$name <- rownames(la2_day2)
la2_day2$countChange <- get.count.change(la2_day2$baseMean, la2_day2$log2FoldChange)
la2_day2 <- la2_day2[order(la2_day2$padj, decreasing = F),]
la2_day2 <- la2_day2[!is.na(la2_day2$padj),]
data.frame(la2_day2)


#>>LA2, DAY4####
se_la2.4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.4$sample]), 
                                 colData = sample.table_la2.4)
ddsSE_2.4 <- DESeqDataSet(se_la2.4, design = ~ con_time)
ddsSE_2.4

#>DESeq run ####
dds_2.4 <- estimateSizeFactors(ddsSE_2.4)
nc <- counts(dds_2.4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.4 <- dds_2.4[filter,]
dds_2.4 <- estimateDispersions(dds_2.4)
dds_2.4 <- nbinomWaldTest(dds_2.4, maxit=500)
resultsNames(dds_2.4)

#4 days vs control
la2_day4 <- results(dds_2.4, contrast = c("con_time", "CNTF_DAY4", "CON_DAY4"))
plotMA(la2_day4)
summary(la2_day4)

la2_day4$name <- rownames(la2_day4)
la2_day4$countChange <- get.count.change(la2_day4$baseMean, la2_day4$log2FoldChange)
la2_day4 <- la2_day4[order(la2_day4$padj, decreasing = F),]
la2_day4 <- la2_day4[!is.na(la2_day4$padj),]
data.frame(la2_day4)




#>>LA5, DAY2####
se_la5.2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la5.2$sample]), 
                                 colData = sample.table_la5.2)
ddsSE_5.2 <- DESeqDataSet(se_la5.2, design = ~ con_time)
ddsSE_5.2

#>DESeq run ####
dds_5.2 <- estimateSizeFactors(ddsSE_5.2)
nc <- counts(dds_5.2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_5.2 <- dds_5.2[filter,]
dds_5.2 <- estimateDispersions(dds_5.2)
dds_5.2 <- nbinomWaldTest(dds_5.2, maxit=500)
resultsNames(dds_5.2)

#30 minutes vs control
la5_min30 <- results(dds_5.2, contrast = c("con_time", "CNTF_MIN30", "CON_DAY2"))
plotMA(la5_min30)
summary(la5_min30)

la5_min30$name <- rownames(la5_min30)
la5_min30$countChange <- get.count.change(la5_min30$baseMean, la5_min30$log2FoldChange)
la5_min30 <- la5_min30[order(la5_min30$padj, decreasing = F),]
la5_min30 <- la5_min30[!is.na(la5_min30$padj),]
data.frame(la5_min30)


#60 minutes vs control
la5_min60 <- results(dds_5.2, contrast = c("con_time", "CNTF_MIN60", "CON_DAY2"))
plotMA(la5_min60)
summary(la5_min60)

la5_min60$name <- rownames(la5_min60)
la5_min60$countChange <- get.count.change(la5_min60$baseMean, la5_min60$log2FoldChange)
la5_min60 <- la5_min60[order(la5_min60$padj, decreasing = F),]
la5_min60 <- la5_min60[!is.na(la5_min60$padj),]
data.frame(la5_min60)

#2 days vs control
la5_day2 <- results(dds_5.2, contrast = c("con_time", "CNTF_DAY2", "CON_DAY2"))
plotMA(la5_day2)
summary(la5_day2)

la5_day2$name <- rownames(la5_day2)
la5_day2$countChange <- get.count.change(la5_day2$baseMean, la5_day2$log2FoldChange)
la5_day2 <- la5_day2[order(la5_day2$padj, decreasing = F),]
la5_day2 <- la5_day2[!is.na(la5_day2$padj),]
data.frame(la5_day2)


#>>LA5, DAY4####
se_la5.4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la5.4$sample]), 
                                 colData = sample.table_la5.4)
ddsSE_5.4 <- DESeqDataSet(se_la5.4, design = ~ con_time)
ddsSE_5.4

#>DESeq run ####
dds_5.4 <- estimateSizeFactors(ddsSE_5.4)
nc <- counts(dds_5.4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_5.4 <- dds_5.4[filter,]
dds_5.4 <- estimateDispersions(dds_5.4)
dds_5.4 <- nbinomWaldTest(dds_5.4, maxit=500)
resultsNames(dds_5.4)

#4 days vs control
la5_day4 <- results(dds_5.4, contrast = c("con_time", "CNTF_DAY4", "CON_DAY4"))
plotMA(la5_day4)
summary(la5_day4)

la5_day4$name <- rownames(la5_day4)
la5_day4$countChange <- get.count.change(la5_day4$baseMean, la5_day4$log2FoldChange)
la5_day4 <- la5_day4[order(la5_day4$padj, decreasing = F),]
la5_day4 <- la5_day4[!is.na(la5_day4$padj),]
data.frame(la5_day4)


#INTERACTION TERM OF CNTF RESPONSE BETWEEN LA2 AND LA5####
sample.table_day2 <- sample.table[sample.table$time != "DAY4",]
sample.table_day2$con_time <- as.factor(as.character(sample.table_day2$con_time))
sample.table_day2$con_time <- relevel(sample.table_day2$con_time, "CON_DAY2")
sample.table_day4 <- sample.table[sample.table$time == "DAY4",]
sample.table_day4$con_time <- as.factor(as.character(sample.table_day4$con_time))
sample.table_day4$con_time <- relevel(sample.table_day4$con_time, "CON_DAY4")

#>>DAY2####
se_la2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_day2$sample]), 
                               colData = sample.table_day2)
ddsSE_2 <- DESeqDataSet(se_la2, design = ~ type + con_time + type:con_time)
ddsSE_2

#>DESeq run ####
dds_2 <- estimateSizeFactors(ddsSE_2)
nc <- counts(dds_2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2 <- dds_2[filter,]
dds_2 <- estimateDispersions(dds_2)
dds_2 <- nbinomWaldTest(dds_2, maxit=500)
resultsNames(dds_2)

#30 minutes vs control
la2vsla5_min30 <- results(dds_2, contrast = list(c("con_time_CNTF_MIN30_vs_CON_DAY2"),
                                                 c("typeLA5.con_timeCNTF_MIN30")))
plotMA(la2vsla5_min30)
summary(la2vsla5_min30)

la2vsla5_min30$name <- rownames(la2vsla5_min30)
la2vsla5_min30$countChange <- get.count.change(la2vsla5_min30$baseMean, la2vsla5_min30$log2FoldChange)
la2vsla5_min30 <- la2vsla5_min30[order(la2vsla5_min30$padj, decreasing = F),]
la2vsla5_min30 <- la2vsla5_min30[!is.na(la2vsla5_min30$padj),]
data.frame(la2vsla5_min30)

#60 minutes vs control
la2vsla5_min60 <- results(dds_2, contrast = list(c("con_time_CNTF_MIN60_vs_CON_DAY2"),
                                                 c("typeLA5.con_timeCNTF_MIN60")))
plotMA(la2vsla5_min60)
summary(la2vsla5_min60)

la2vsla5_min60$name <- rownames(la2vsla5_min60)
la2vsla5_min60$countChange <- get.count.change(la2vsla5_min60$baseMean, la2vsla5_min60$log2FoldChange)
la2vsla5_min60 <- la2vsla5_min60[order(la2vsla5_min60$padj, decreasing = F),]
la2vsla5_min60 <- la2vsla5_min60[!is.na(la2vsla5_min60$padj),]
data.frame(la2vsla5_min60)


#2 days vs control
la2vsla5_day2 <- results(dds_2, contrast = list(c("con_time_CNTF_DAY2_vs_CON_DAY2"),
                                                c("typeLA5.con_timeCNTF_DAY2")))
plotMA(la2vsla5_day2)
summary(la2vsla5_day2)

la2vsla5_day2$name <- rownames(la2vsla5_day2)
la2vsla5_day2$countChange <- get.count.change(la2vsla5_day2$baseMean, la2vsla5_day2$log2FoldChange)
la2vsla5_day2 <- la2vsla5_day2[order(la2vsla5_day2$padj, decreasing = F),]
la2vsla5_day2 <- la2vsla5_day2[!is.na(la2vsla5_day2$padj),]
data.frame(la2vsla5_day2)


#>>DAY4####
se_la4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_day4$sample]), 
                               colData = sample.table_day4)
ddsSE_4 <- DESeqDataSet(se_la4, design = ~ type + con_time + type:con_time)
ddsSE_4

#>DESeq run ####
dds_4 <- estimateSizeFactors(ddsSE_4)
nc <- counts(dds_4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_4 <- dds_4[filter,]
dds_4 <- estimateDispersions(dds_4)
dds_4 <- nbinomWaldTest(dds_4, maxit=500)
resultsNames(dds_4)

#4 days vs control
la2vsla5_day4 <- results(dds_4, contrast = list(c("con_time_CNTF_DAY4_vs_CON_DAY4"),
                                                c("typeLA5.con_timeCNTF_DAY4")))
plotMA(la2vsla5_day4)
summary(la2vsla5_day4)

la2vsla5_day4$name <- rownames(la2vsla5_day4)
la2vsla5_day4$countChange <- get.count.change(la2vsla5_day4$baseMean, la2vsla5_day4$log2FoldChange)
la2vsla5_day4 <- la2vsla5_day4[order(la2vsla5_day4$padj, decreasing = F),]
la2vsla5_day4 <- la2vsla5_day4[!is.na(la2vsla5_day4$padj),]
data.frame(la2vsla5_day4)

#save####
save(la2_day2, la2_day4, la2_min30, la2_min60, 
     la5_day2, la5_day4, la5_min30, la5_min60, 
     la2vsla5_day2, la2vsla5_day4, la2vsla5_min30, la2vsla5_min60,
     file = "working_data/mirna_seq/LA2_LA5_DESeq_results_new_seq.RData")

#EXPERIMENT LA2 AND LA2_LA5 + BATCHES####
#new seq
mirexp.folder <- "./working_data/mirna_seq/mirexpress_alignment/LA2_LA5/results"
mirexp.samples <- sort(list.dirs(mirexp.folder, full.names = F, recursive = F))
mirexp.dirs <- sort(list.dirs(mirexp.folder, recursive = F))

n <- length(mirexp.dirs)
mirexp <- vector(mode = "list", length = n)
names(mirexp) <- mirexp.samples
for(i in 1:n){
  temp <- read.table(paste0(mirexp.dirs[i], "/miRNA_expression.txt"), sep = "\t")
  colnames(temp)[2] <- "count"
  temp$anc <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 1))
  temp$mat <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 2))
  temp$V1 <- NULL
  head(temp)
  temp <- temp[order(temp$count, decreasing = T),]
  temp <- temp[which(!duplicated(temp$mat)),]
  mirexp[[i]] <- temp
}

mirexp_new <- mirexp

#old seq
mirexp.folder <- "./working_data/mirna_seq/mirexpress_alignment/LA2"
mirexp.files <- list.files(mirexp.folder, full.names = T)
mirexp.files <- mirexp.files[grep("R1.*expression.txt", mirexp.files)]

n <- length(mirexp.files)
mirexp <- vector(mode = "list", length = n)
for(i in 1:n){
  temp <- read.table(mirexp.files[i], sep = "\t")
  colnames(temp)[2] <- "count"
  temp$anc <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 1))
  temp$mat <- unlist(lapply(strsplit(as.character(temp$V1), ","), "[", 2))
  temp$V1 <- NULL
  head(temp)
  temp <- temp[order(temp$count, decreasing = T),]
  temp <- temp[which(!duplicated(temp$mat)),]
  mirexp[[i]] <- temp
}

mirexp_old <- mirexp
names(mirexp_old) <- seq(18)

mirexp <- c(mirexp_new, mirexp_old)

all_mirs <- sort(unique(unlist(lapply(mirexp, function(x) x$mat))))
counts <- data.frame(name = all_mirs)
n <- length(mirexp)
for(i in 1:n) {
  idx <- match(counts$name, mirexp[[i]]$mat)
  counts[, (i+1)] <- mirexp[[i]]$count[idx]
}

rownames(counts) <- counts$name
counts$name <- NULL
colnames(counts) <- names(mirexp)
counts[is.na(counts)] <- 0
head(counts)

colnames(counts) <- paste0("sample", formatC(seq(1, 24), width = 2, format = "d", flag = "0"))
write.csv(counts[,1:18], file = "out/mirexpress2_counts.csv", row.names = T)

## sample table ####
sample.table <- data.frame(sample = colnames(counts))
rownames(sample.table) <- sample.table$sample
sample.table$sample[49:nrow(sample.table)] <- paste(49:nrow(sample.table), "-",
                                    rep("LA2", 18),
                                    c(rep("CO", 3), rep("CN", 3), rep("CO", 3), rep("CN", 3), rep("CN", 6)),
                                    c(rep(2, 6), rep(4, 6), rep(30, 3), rep(60, 3)), "-",
                                    rep(5:7),
                                    sep = "")
sample.table$sample.name <- unlist(lapply(strsplit((sample.table$sample), "-"), "[", 2))
sample.table$condition <- as.factor(grepl("CN", sample.table$sample.name))
levels(sample.table$condition) <- c("CON", "CNTF")

sample.table$batch <- as.factor(c(rep(1, 48), rep(2, 18)))

sample.table$time <- sample.table$sample.name
sample.table$time <- gsub("LA2CO", "", sample.table$time); sample.table$time <- gsub("LA5CO", "", sample.table$time); 
sample.table$time <- gsub("LA2CN", "", sample.table$time); sample.table$time <- gsub("LA5CN", "", sample.table$time); 
sample.table$time <- as.factor(sample.table$time)
levels(sample.table$time) <- c("DAY2", "MIN30", "DAY4", "MIN60")

sample.table$type <- as.factor(substr(sample.table$sample.name, 1, 3))

sample.table$con_time_type <- as.factor(apply(sample.table[, c("condition", "time", "type")], 1, paste, collapse = "_"))
sample.table$con_time <- as.factor(apply(sample.table[, c("condition", "time")], 1, paste, collapse = "_"))

dim(counts)
dim(sample.table)
all(colnames(counts) == rownames(sample.table))
rownames(sample.table) <- colnames(counts) <- sample.table$sample

#split sample tables per cell type and day for simple DE####
sample.table_la2 <- sample.table[sample.table$type == "LA2",]
sample.table_la2.2 <- sample.table_la2[sample.table_la2$time != "DAY4",]
sample.table_la2.2$con_time <- as.factor(as.character(sample.table_la2.2$con_time))
sample.table_la2.2$con_time <- relevel(sample.table_la2.2$con_time, "CON_DAY2")
sample.table_la2.4 <- sample.table_la2[sample.table_la2$time == "DAY4",]
sample.table_la2.4$con_time <- as.factor(as.character(sample.table_la2.4$con_time))
sample.table_la2.4$con_time <- relevel(sample.table_la2.4$con_time, "CON_DAY4")

sample.table_la5 <- sample.table[sample.table$type == "LA5",]
sample.table_la5.2 <- sample.table_la5[sample.table_la5$time != "DAY4",]
sample.table_la5.2$con_time <- as.factor(as.character(sample.table_la5.2$con_time))
sample.table_la5.2$con_time <- relevel(sample.table_la5.2$con_time, "CON_DAY2")
sample.table_la5.4 <- sample.table_la5[sample.table_la5$time == "DAY4",]
sample.table_la5.4$con_time <- as.factor(as.character(sample.table_la5.4$con_time))
sample.table_la5.4$con_time <- relevel(sample.table_la5.4$con_time, "CON_DAY4")


## DESeq SINGLE SETS DE ####
#>>LA2, DAY2####
se_la2.2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.2$sample]), 
                                 colData = sample.table_la2.2)
ddsSE_2.2 <- DESeqDataSet(se_la2.2, design = ~ batch + con_time)
ddsSE_2.2

#>DESeq run ####
dds_2.2 <- estimateSizeFactors(ddsSE_2.2)
nc <- counts(dds_2.2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.2 <- dds_2.2[filter,]
dds_2.2 <- estimateDispersions(dds_2.2)
dds_2.2 <- nbinomWaldTest(dds_2.2, maxit=500)
resultsNames(dds_2.2)

#30 minutes vs control
la2_min30 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN30", "CON_DAY2"))
plotMA(la2_min30)
summary(la2_min30)

la2_min30$name <- rownames(la2_min30)
la2_min30$countChange <- get.count.change(la2_min30$baseMean, la2_min30$log2FoldChange)
la2_min30 <- la2_min30[order(la2_min30$padj, decreasing = F),]
la2_min30 <- la2_min30[!is.na(la2_min30$padj),]
data.frame(la2_min30)


#60 minutes vs control
la2_min60 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN60", "CON_DAY2"))
plotMA(la2_min60)
summary(la2_min60)

la2_min60$name <- rownames(la2_min60)
la2_min60$countChange <- get.count.change(la2_min60$baseMean, la2_min60$log2FoldChange)
la2_min60 <- la2_min60[order(la2_min60$padj, decreasing = F),]
la2_min60 <- la2_min60[!is.na(la2_min60$padj),]
data.frame(la2_min60)

#2 days vs control
la2_day2 <- results(dds_2.2, contrast = c("con_time", "CNTF_DAY2", "CON_DAY2"))
plotMA(la2_day2)
summary(la2_day2)

la2_day2$name <- rownames(la2_day2)
la2_day2$countChange <- get.count.change(la2_day2$baseMean, la2_day2$log2FoldChange)
la2_day2 <- la2_day2[order(la2_day2$padj, decreasing = F),]
la2_day2 <- la2_day2[!is.na(la2_day2$padj),]
data.frame(la2_day2)


#>>LA2, DAY4####
se_la2.4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.4$sample]), 
                                 colData = sample.table_la2.4)
ddsSE_2.4 <- DESeqDataSet(se_la2.4, design = ~ batch + con_time)
ddsSE_2.4

#>DESeq run ####
dds_2.4 <- estimateSizeFactors(ddsSE_2.4)
nc <- counts(dds_2.4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.4 <- dds_2.4[filter,]
dds_2.4 <- estimateDispersions(dds_2.4)
dds_2.4 <- nbinomWaldTest(dds_2.4, maxit=500)
resultsNames(dds_2.4)

#4 days vs control
la2_day4 <- results(dds_2.4, contrast = c("con_time", "CNTF_DAY4", "CON_DAY4"))
plotMA(la2_day4)
summary(la2_day4)

la2_day4$name <- rownames(la2_day4)
la2_day4$countChange <- get.count.change(la2_day4$baseMean, la2_day4$log2FoldChange)
la2_day4 <- la2_day4[order(la2_day4$padj, decreasing = F),]
la2_day4 <- la2_day4[!is.na(la2_day4$padj),]
data.frame(la2_day4)

#>LA5, should be unchanged####
#>>LA5, DAY2####
se_la5.2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la5.2$sample]), 
                                 colData = sample.table_la5.2)
ddsSE_5.2 <- DESeqDataSet(se_la5.2, design = ~ batch + con_time)
ddsSE_5.2

#>DESeq run ####
dds_5.2 <- estimateSizeFactors(ddsSE_5.2)
nc <- counts(dds_5.2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_5.2 <- dds_5.2[filter,]
dds_5.2 <- estimateDispersions(dds_5.2)
dds_5.2 <- nbinomWaldTest(dds_5.2, maxit=500)
resultsNames(dds_5.2)

#30 minutes vs control
la5_min30 <- results(dds_5.2, contrast = c("con_time", "CNTF_MIN30", "CON_DAY2"))
plotMA(la5_min30)
summary(la5_min30)

la5_min30$name <- rownames(la5_min30)
la5_min30$countChange <- get.count.change(la5_min30$baseMean, la5_min30$log2FoldChange)
la5_min30 <- la5_min30[order(la5_min30$padj, decreasing = F),]
la5_min30 <- la5_min30[!is.na(la5_min30$padj),]
data.frame(la5_min30)


#60 minutes vs control
la5_min60 <- results(dds_5.2, contrast = c("con_time", "CNTF_MIN60", "CON_DAY2"))
plotMA(la5_min60)
summary(la5_min60)

la5_min60$name <- rownames(la5_min60)
la5_min60$countChange <- get.count.change(la5_min60$baseMean, la5_min60$log2FoldChange)
la5_min60 <- la5_min60[order(la5_min60$padj, decreasing = F),]
la5_min60 <- la5_min60[!is.na(la5_min60$padj),]
data.frame(la5_min60)

#2 days vs control
la5_day2 <- results(dds_5.2, contrast = c("con_time", "CNTF_DAY2", "CON_DAY2"))
plotMA(la5_day2)
summary(la5_day2)

la5_day2$name <- rownames(la5_day2)
la5_day2$countChange <- get.count.change(la5_day2$baseMean, la5_day2$log2FoldChange)
la5_day2 <- la5_day2[order(la5_day2$padj, decreasing = F),]
la5_day2 <- la5_day2[!is.na(la5_day2$padj),]
data.frame(la5_day2)


#>>LA5, DAY4####
se_la5.4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la5.4$sample]), 
                                 colData = sample.table_la5.4)
ddsSE_5.4 <- DESeqDataSet(se_la5.4, design = ~ con_time)
colnames(ddsSE_5.4)

#>DESeq run ####
dds_5.4 <- estimateSizeFactors(ddsSE_5.4)
nc <- counts(dds_5.4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_5.4 <- dds_5.4[filter,]
dds_5.4 <- estimateDispersions(dds_5.4)
dds_5.4 <- nbinomWaldTest(dds_5.4, maxit=500)
resultsNames(dds_5.4)

#4 days vs control
la5_day4 <- results(dds_5.4, contrast = c("con_time", "CNTF_DAY4", "CON_DAY4"))
plotMA(la5_day4)
summary(la5_day4)

la5_day4$name <- rownames(la5_day4)
la5_day4$countChange <- get.count.change(la5_day4$baseMean, la5_day4$log2FoldChange)
la5_day4 <- la5_day4[order(la5_day4$padj, decreasing = F),]
la5_day4 <- la5_day4[!is.na(la5_day4$padj),]
data.frame(la5_day4)


#INTERACTION TERM OF CNTF RESPONSE BETWEEN LA2 AND LA5####
sample.table_day2 <- sample.table[sample.table$time != "DAY4",]
sample.table_day2$con_time <- as.factor(as.character(sample.table_day2$con_time))
sample.table_day2$con_time <- relevel(sample.table_day2$con_time, "CON_DAY2")
sample.table_day4 <- sample.table[sample.table$time == "DAY4",]
sample.table_day4$con_time <- as.factor(as.character(sample.table_day4$con_time))
sample.table_day4$con_time <- relevel(sample.table_day4$con_time, "CON_DAY4")

#>>DAY2####
se_la2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_day2$sample]), 
                               colData = sample.table_day2)
ddsSE_2 <- DESeqDataSet(se_la2, design = ~ batch + type + con_time + type:con_time)
ddsSE_2

#>DESeq run ####
dds_2 <- estimateSizeFactors(ddsSE_2)
nc <- counts(dds_2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2 <- dds_2[filter,]
dds_2 <- estimateDispersions(dds_2)
dds_2 <- nbinomWaldTest(dds_2, maxit=500)
resultsNames(dds_2)

#30 minutes vs control
la2vsla5_min30 <- results(dds_2, contrast = list(c("con_time_CNTF_MIN30_vs_CON_DAY2"),
                                                 c("typeLA5.con_timeCNTF_MIN30")))
plotMA(la2vsla5_min30)
summary(la2vsla5_min30)

la2vsla5_min30$name <- rownames(la2vsla5_min30)
la2vsla5_min30$countChange <- get.count.change(la2vsla5_min30$baseMean, la2vsla5_min30$log2FoldChange)
la2vsla5_min30 <- la2vsla5_min30[order(la2vsla5_min30$padj, decreasing = F),]
la2vsla5_min30 <- la2vsla5_min30[!is.na(la2vsla5_min30$padj),]
data.frame(la2vsla5_min30)

#60 minutes vs control
la2vsla5_min60 <- results(dds_2, contrast = list(c("con_time_CNTF_MIN60_vs_CON_DAY2"),
                                                 c("typeLA5.con_timeCNTF_MIN60")))
plotMA(la2vsla5_min60)
summary(la2vsla5_min60)

la2vsla5_min60$name <- rownames(la2vsla5_min60)
la2vsla5_min60$countChange <- get.count.change(la2vsla5_min60$baseMean, la2vsla5_min60$log2FoldChange)
la2vsla5_min60 <- la2vsla5_min60[order(la2vsla5_min60$padj, decreasing = F),]
la2vsla5_min60 <- la2vsla5_min60[!is.na(la2vsla5_min60$padj),]
data.frame(la2vsla5_min60)


#2 days vs control
la2vsla5_day2 <- results(dds_2, contrast = list(c("con_time_CNTF_DAY2_vs_CON_DAY2"),
                                                c("typeLA5.con_timeCNTF_DAY2")))
plotMA(la2vsla5_day2)
summary(la2vsla5_day2)

la2vsla5_day2$name <- rownames(la2vsla5_day2)
la2vsla5_day2$countChange <- get.count.change(la2vsla5_day2$baseMean, la2vsla5_day2$log2FoldChange)
la2vsla5_day2 <- la2vsla5_day2[order(la2vsla5_day2$padj, decreasing = F),]
la2vsla5_day2 <- la2vsla5_day2[!is.na(la2vsla5_day2$padj),]
data.frame(la2vsla5_day2)


#>>DAY4####
se_la4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_day4$sample]), 
                               colData = sample.table_day4)
ddsSE_4 <- DESeqDataSet(se_la4, design = ~ batch + type + con_time + type:con_time)
ddsSE_4

#>DESeq run ####
dds_4 <- estimateSizeFactors(ddsSE_4)
nc <- counts(dds_4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_4 <- dds_4[filter,]
dds_4 <- estimateDispersions(dds_4)
dds_4 <- nbinomWaldTest(dds_4, maxit=500)
resultsNames(dds_4)

#4 days vs control
la2vsla5_day4 <- results(dds_4, contrast = list(c("con_time_CNTF_DAY4_vs_CON_DAY4"),
                                                c("typeLA5.con_timeCNTF_DAY4")))
plotMA(la2vsla5_day4)
summary(la2vsla5_day4)

la2vsla5_day4$name <- rownames(la2vsla5_day4)
la2vsla5_day4$countChange <- get.count.change(la2vsla5_day4$baseMean, la2vsla5_day4$log2FoldChange)
la2vsla5_day4 <- la2vsla5_day4[order(la2vsla5_day4$padj, decreasing = F),]
la2vsla5_day4 <- la2vsla5_day4[!is.na(la2vsla5_day4$padj),]
data.frame(la2vsla5_day4)

#save####
save(la2_day2, la2_day4, la2_min30, la2_min60, 
     la5_day2, la5_day4, la5_min30, la5_min60, 
     la2vsla5_day2, la2vsla5_day4, la2vsla5_min30, la2vsla5_min60,
     file = "working_data/mirna_seq/LA2_LA5_DESeq_results_all_seq_batches.RData")

#EXPERIMENT ONLY LA2####
sample.table_la2 <- sample.table[49:66,]
#split sample tables per cell type and day for simple DE####
sample.table_la2.2 <- sample.table_la2[sample.table_la2$time != "DAY4",]
sample.table_la2.2$con_time <- as.factor(as.character(sample.table_la2.2$con_time))
sample.table_la2.2$con_time <- relevel(sample.table_la2.2$con_time, "CON_DAY2")
sample.table_la2.4 <- sample.table_la2[sample.table_la2$time == "DAY4",]
sample.table_la2.4$con_time <- as.factor(as.character(sample.table_la2.4$con_time))
sample.table_la2.4$con_time <- relevel(sample.table_la2.4$con_time, "CON_DAY4")


## DESeq SINGLE SETS DE ####
#>>LA2, DAY2####
se_la2.2 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.2$sample]), 
                                 colData = sample.table_la2.2)
ddsSE_2.2 <- DESeqDataSet(se_la2.2, design = ~ con_time)
ddsSE_2.2

#>DESeq run ####
dds_2.2 <- estimateSizeFactors(ddsSE_2.2)
nc <- counts(dds_2.2, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.2 <- dds_2.2[filter,]
dds_2.2 <- estimateDispersions(dds_2.2)
dds_2.2 <- nbinomWaldTest(dds_2.2, maxit=500)
resultsNames(dds_2.2)

#30 minutes vs control
la2_min30 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN30", "CON_DAY2"))
plotMA(la2_min30)
summary(la2_min30)

la2_min30$name <- rownames(la2_min30)
la2_min30$countChange <- get.count.change(la2_min30$baseMean, la2_min30$log2FoldChange)
la2_min30 <- la2_min30[order(la2_min30$padj, decreasing = F),]
la2_min30 <- la2_min30[!is.na(la2_min30$padj),]
data.frame(la2_min30)


#60 minutes vs control
la2_min60 <- results(dds_2.2, contrast = c("con_time", "CNTF_MIN60", "CON_DAY2"))
plotMA(la2_min60)
summary(la2_min60)

la2_min60$name <- rownames(la2_min60)
la2_min60$countChange <- get.count.change(la2_min60$baseMean, la2_min60$log2FoldChange)
la2_min60 <- la2_min60[order(la2_min60$padj, decreasing = F),]
la2_min60 <- la2_min60[!is.na(la2_min60$padj),]
data.frame(la2_min60)

#2 days vs control
la2_day2 <- results(dds_2.2, contrast = c("con_time", "CNTF_DAY2", "CON_DAY2"))
plotMA(la2_day2)
summary(la2_day2)

la2_day2$name <- rownames(la2_day2)
la2_day2$countChange <- get.count.change(la2_day2$baseMean, la2_day2$log2FoldChange)
la2_day2 <- la2_day2[order(la2_day2$padj, decreasing = F),]
la2_day2 <- la2_day2[!is.na(la2_day2$padj),]
data.frame(la2_day2)


#>>LA2, DAY4####
se_la2.4 <- SummarizedExperiment(as.matrix(counts[, colnames(counts) %in% sample.table_la2.4$sample]), 
                                 colData = sample.table_la2.4)
ddsSE_2.4 <- DESeqDataSet(se_la2.4, design = ~ con_time)
ddsSE_2.4

#>DESeq run ####
dds_2.4 <- estimateSizeFactors(ddsSE_2.4)
nc <- counts(dds_2.4, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
mean(filter)
dds_2.4 <- dds_2.4[filter,]
dds_2.4 <- estimateDispersions(dds_2.4)
dds_2.4 <- nbinomWaldTest(dds_2.4, maxit=500)
resultsNames(dds_2.4)

#4 days vs control
la2_day4 <- results(dds_2.4, contrast = c("con_time", "CNTF_DAY4", "CON_DAY4"))
plotMA(la2_day4)
summary(la2_day4)

la2_day4$name <- rownames(la2_day4)
la2_day4$countChange <- get.count.change(la2_day4$baseMean, la2_day4$log2FoldChange)
la2_day4 <- la2_day4[order(la2_day4$padj, decreasing = F),]
la2_day4 <- la2_day4[!is.na(la2_day4$padj),]
data.frame(la2_day4)

#save####
save(la2_day2, la2_day4, la2_min30, la2_min60, 
     file = "working_data/mirna_seq/LA2_LA5_DESeq_results_old_seq_validation.RData")
