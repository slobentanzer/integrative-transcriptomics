#NULL DISTRIBUTION FOR PERMUTATION TEST OF SEX SPECIFIC META ANALYSES, 10 000 ITERATIONS####

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
#MALE####
scz_bd_meta_dge_m <- readRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge_m.RData")
genes_m = rownames(scz_bd_meta_dge_m[[1]]$exprs)
compute_null_m <- function(it) {
  # message(it)
  seed <- round(runif(1)*1e7, 0)
  allmeta = matrix(NA, nrow = length(genes_m), length(scz_bd_meta_dge_m))
  colnames(allmeta) = c("SCZ", "BD")
  allmeta =  as.data.frame(allmeta)
  for (i in 1:length(scz_bd_meta_dge_m)) {
    # print(i)
    tt = matrix(NA, nrow = length(genes_m), ncol = 3)
    subj = unique(as.character(scz_bd_meta_dge_m[[i]]$meta$Subject))
    subj_group = data.frame(Subject = subj, Group = scz_bd_meta_dge_m[[i]]$meta$Group[match(subj, scz_bd_meta_dge_m[[i]]$meta$Subject)])
    subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
    scz_bd_meta_dge_m[[i]]$meta$Group = subj_group$Group[match(scz_bd_meta_dge_m[[i]]$meta$Subject, subj_group$Subject)]
    
    allmeta[, i] = apply(scz_bd_meta_dge_m[[i]]$exprs, 1, lmer_apply, scz_bd_meta_dge_m[[i]]$meta)
    
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
    file = paste0("working_data/01_unbiased/null_dist_m/", seed, ".txt"),
    row.names = F,
    col.names = F
  )
}


#FEMALE####
scz_bd_meta_dge_f <- readRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge_f.RData")
genes_f = rownames(scz_bd_meta_dge_f[[1]]$exprs)

compute_null_f <- function(it) {
  # message(it)
  seed <- round(runif(1)*1e7, 0)
  allmeta = matrix(NA, nrow = length(genes_f), length(scz_bd_meta_dge_f))
  colnames(allmeta) = c("SCZ", "BD")
  allmeta =  as.data.frame(allmeta)
  for (i in 1:length(scz_bd_meta_dge_f)) {
    # print(i)
    tt = matrix(NA, nrow = length(genes_f), ncol = 3)
    subj = unique(as.character(scz_bd_meta_dge_f[[i]]$meta$Subject))
    subj_group = data.frame(Subject = subj, Group = scz_bd_meta_dge_f[[i]]$meta$Group[match(subj, scz_bd_meta_dge_f[[i]]$meta$Subject)])
    subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
    scz_bd_meta_dge_f[[i]]$meta$Group = subj_group$Group[match(scz_bd_meta_dge_f[[i]]$meta$Subject, subj_group$Subject)]
    
    allmeta[, i] = apply(scz_bd_meta_dge_f[[i]]$exprs, 1, lmer_apply, scz_bd_meta_dge_f[[i]]$meta)
    
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
    file = paste0("working_data/01_unbiased/null_dist_f/", seed, ".txt"),
    row.names = F,
    col.names = F
  )
}

#PARALLEL APPLICATION####
#male
files=paste0("working_data/01_unbiased/null_dist_m/", dir(path = "working_data/01_unbiased/null_dist_m/", pattern="*.txt"))
seq <- seq((10000 - length(files))/6) #6 processes
# bplapply(seq, compute_null_m, BPPARAM=MulticoreParam(6, tasks=length(seq), progressbar = T))
# system.time(compute_null_m(1))
for(i in seq){
  print(i)
  compute_null_m(i)
}

#female
files=paste0("working_data/unbiased/null_dist_f/", dir(path = "working_data/unbiased/null_dist_f/", pattern="*.txt"))
seq <- seq((10000 - length(files))/6)
# bplapply(seq, compute_null_f, BPPARAM=MulticoreParam(6, tasks=length(seq), progressbar = T))
for(i in seq){
  print(i)
  compute_null_f(i)
}


#READ####
svg(filename = "./img/unbiased_null_dist.svg", width = 10, height = 3)
par(mfrow = c(1,3))
#independent
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
hist(null$null, 50, main="Null Distribution\n10,000 permutations", xlab="Spearman's Rho", breaks = 20)

#male
files=paste0("working_data/unbiased/null_dist_m/", dir(path = "working_data/unbiased/null_dist_m/", pattern="*.txt"))
null = vector(mode="numeric")

for(i in 1:length(files)) {
  d = read.table(files[[i]])
  null = c(null, d$V1)
}
print(length(null))
print(quantile(null,probs = c(0,.025, .5, .975, 1)))

null = sort(null)
null = cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
hist(null$null, 50, main="Null Distribution, males\n10,000 permutations", xlab="Spearman's Rho", breaks = 20)

bd_meta_dge_m <- readRDS(file = "out/bd_meta_dge_m.rds")
scz_meta_dge_m <- readRDS(file = "out/scz_meta_dge_m.rds")
R <- cor.test(scz_meta_dge_m$beta, bd_meta_dge_m$beta, method = "spear")
rho <- cor(scz_meta_dge_m$beta, bd_meta_dge_m$beta, method = "spear", use = "complete")
sem <- (tanh(atanh(rho + 1.96/sqrt(nrow(scz_meta_dge_m)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(scz_meta_dge_m)-3))))/3.92
p.bootstrap <- null$prob[findInterval(rho, null$null)]

R; rho; sem; p.bootstrap

#female
files=paste0("working_data/unbiased/null_dist_f/", dir(path = "working_data/unbiased/null_dist_f/", pattern="*.txt"))
null = vector(mode="numeric")

for(i in 1:length(files)) {
  d = read.table(files[[i]])
  null = c(null, d$V1)
}
print(length(null))
print(quantile(null,probs = c(0,.025, .5, .975, 1)))

null = sort(null)
null = cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
hist(null$null, 50, main="Null Distribution, females\n10,000 permutations", xlab="Spearman's Rho", breaks = 20)
dev.off()

bd_meta_dge_f <- readRDS(file = "out/bd_meta_dge_f.rds")
scz_meta_dge_f <- readRDS(file = "out/scz_meta_dge_f.rds")
R <- cor.test(scz_meta_dge_f$beta, bd_meta_dge_f$beta, method = "spear")
rho <- cor(scz_meta_dge_f$beta, bd_meta_dge_f$beta, method = "spear", use = "complete")
sem <- (tanh(atanh(rho + 1.96/sqrt(nrow(scz_meta_dge_f)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(scz_meta_dge_f)-3))))/3.92
p.bootstrap <- null$prob[findInterval(rho, null$null)]

R; rho; sem; p.bootstrap
