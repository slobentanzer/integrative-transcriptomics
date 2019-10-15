#NULL DISTRIBUTION FOR PERMUTATION TEST OF FULL META ANALYSIS, 10 000 ITERATIONS####

rm(list=ls())
library(BiocParallel)
registered()
library(lme4)

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

lmer_apply=function(x, meta) {
  if(length(unique(meta$Subject)) < nrow(meta)) {
    # message("A")
    return(summary(lmer(x ~ Group + Study + (1 | Subject),data=meta))$coefficients[2,1])
  } else {
    # message("B")
    return(summary(lm(x ~ Group + Study,data=meta))$coefficients[2,1])
  }
}

#LOAD####
scz_bd_meta_dge <- readRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge.RData")
genes = rownames(scz_bd_meta_dge[[1]]$exprs)
compute_null <- function(it) {
  # message(it)
  seed <- round(runif(1)*1e7, 0)
  allmeta = matrix(NA, nrow = length(genes), length(scz_bd_meta_dge))
  colnames(allmeta) = c("SCZ", "BD")
  allmeta =  as.data.frame(allmeta)
  for (i in 1:length(scz_bd_meta_dge)) {
    # print(i)
    tt = matrix(NA, nrow = length(genes), ncol = 3)
    subj = unique(as.character(scz_bd_meta_dge[[i]]$meta$Subject))
    subj_group = data.frame(Subject = subj, Group = scz_bd_meta_dge[[i]]$meta$Group[match(subj, scz_bd_meta_dge[[i]]$meta$Subject)])
    subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
    scz_bd_meta_dge[[i]]$meta$Group = subj_group$Group[match(scz_bd_meta_dge[[i]]$meta$Subject, subj_group$Subject)]
    
    allmeta[, i] = apply(scz_bd_meta_dge[[i]]$exprs, 1, lmer_apply, scz_bd_meta_dge[[i]]$meta)
    
  }
  
  
  cor_vec = vector(mode = "numeric")
  comparisons = t(combn(seq(1, ncol(allmeta)), 2))
  
  for (i in 1:nrow(comparisons)) {
    r = cor(allmeta[, comparisons[i, 1]], allmeta[, comparisons[i, 2]], method = "spearman", use =
              "pairwise.complete.obs")
    cor_vec = c(cor_vec, r)
  }
  
  write.table(
    cor_vec,
    file = paste0("working_data/01_unbiased/null_dist/", seed, ".txt"),
    row.names = F,
    col.names = F
  )
}


#PARALLEL APPLICATION####
files=paste0("working_data/01_unbiased/null_dist/", dir(path = "working_data/01_unbiased/null_dist/", pattern="*.txt"))
seq <- seq((10000 - length(files))/6) #6 processes
# bplapply(seq, compute_null, BPPARAM=MulticoreParam(6, tasks=length(seq), progressbar = T))
# system.time(compute_null(1))
for(i in seq){
  print(i)
  compute_null(i)
}


#READ####
f=paste0("working_data/unbiased/null_dist/", dir(path = "working_data/unbiased/null_dist/", pattern="*.txt"))
null = vector(mode="numeric")

for(i in 1:length(f)) {
  d = read.table(f[[i]])
  null = c(null, d$V1)
}
print(length(null))
print(quantile(null,probs = c(0,.025, .5, .975, 1)))

null <-  sort(null)
null <-  cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
hist(null$null, 50, main="Null Distribution\n10,000 permutations", xlab="Spearman's Rho")

bd_meta_dge <- readRDS(file = "out/bd_meta_dge.rds")
scz_meta_dge <- readRDS(file = "out/scz_meta_dge.rds")
R <- cor.test(scz_meta_dge$beta, bd_meta_dge$beta, method = "spear")
rho <- cor(scz_meta_dge$beta, bd_meta_dge$beta, method = "spear", use = "complete")
sem <- (tanh(atanh(rho + 1.96/sqrt(nrow(scz_meta_dge)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(scz_meta_dge)-3))))/3.92
p.bootstrap <- null$prob[findInterval(rho, null$null)]

R; rho; sem; p.bootstrap
