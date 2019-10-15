#SEX SPECIFIC META-ANALYSES OF BD AND SCZ SAMPLES, EXTENDING GANDAL ET AL####

rm(list=ls())
options(stringsAsFactors=F)

library(nlme)

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#SCZ####
files = dir("./working_data/01_unbiased/", pattern="_normalized")
files = files[grep("sex",files)]
n = length(files)
SCZ_exprs_list = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/01_unbiased/", files[[i]], sep=""))
  SCZ_exprs_list[[i]]$exprs = exprs 
  SCZ_exprs_list[[i]]$meta = meta
  SCZ_exprs_list[[i]]$annot = annot
  rm(exprs,meta,annot)
}

genes = intersect(rownames(SCZ_exprs_list[[1]]$exprs), rownames(SCZ_exprs_list[[2]]$exprs))
for(i in 3:n) genes = intersect(genes, rownames(SCZ_exprs_list[[i]]$exprs))
for(i in 1:n) SCZ_exprs_list[[i]]$exprs = SCZ_exprs_list[[i]]$exprs[match(genes,rownames(SCZ_exprs_list[[i]]$exprs)),]

SCZ_exprs_m = data.frame(row.names = genes)
SCZ_meta_m = data.frame(Study=NA, Subject=NA, Group=NA)
SCZ_exprs_f = data.frame(row.names = genes)
SCZ_meta_f = data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  idx <- which(SCZ_exprs_list[[i]]$meta[,"Sex"] == "M")

  SCZ_exprs_m = cbind(SCZ_exprs_m, SCZ_exprs_list[[i]]$exprs[,idx])
  SCZ_meta_m = rbind(SCZ_meta_m, SCZ_exprs_list[[i]]$meta[idx,c("Study", "Subject","Group")])
  
  SCZ_exprs_f = cbind(SCZ_exprs_f, SCZ_exprs_list[[i]]$exprs[,-idx])
  SCZ_meta_f = rbind(SCZ_meta_f, SCZ_exprs_list[[i]]$meta[-idx,c("Study", "Subject","Group")])
  
}
SCZ_meta_m = SCZ_meta_m[-1,]
to_keep = SCZ_meta_m$Group %in% c("CTL", "SCZ")
SCZ_exprs_m = SCZ_exprs_m[,to_keep]; SCZ_meta_m = SCZ_meta_m[to_keep,]
SCZ_meta_m$Group = factor(SCZ_meta_m$Group)
SCZ_meta_m$Study = as.factor(SCZ_meta_m$Study)
SCZ_meta_m$Subject = as.factor(SCZ_meta_m$Subject)

SCZ_meta_f = SCZ_meta_f[-1,]
to_keep = SCZ_meta_f$Group %in% c("CTL", "SCZ")
SCZ_exprs_f = SCZ_exprs_f[,to_keep]; SCZ_meta_f = SCZ_meta_f[to_keep,]
SCZ_meta_f$Group = factor(SCZ_meta_f$Group)
SCZ_meta_f$Study = as.factor(SCZ_meta_f$Study)
SCZ_meta_f$Subject = as.factor(SCZ_meta_f$Subject)

scz_meta_dge_m = matrix(NA, nrow=nrow(SCZ_exprs_m), ncol=3)
scz_meta_dge_f = matrix(NA, nrow=nrow(SCZ_exprs_f), ncol=3)
for(i in 1:nrow(SCZ_exprs_m)) {
  if(i%%100==0) print(i)
  expr_m = as.numeric(SCZ_exprs_m[i,])
  expr_f = as.numeric(SCZ_exprs_f[i,])
  tryCatch({
    scz_meta_dge_m[i,] = summary(lme(expr_m~ Group + Study,data = SCZ_meta_m, random=~1|Subject))$tTable[2,c(1,2,5)]
    scz_meta_dge_f[i,] = summary(lme(expr_f~ Group + Study,data = SCZ_meta_f, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e) print(e))
}

scz_meta_dge_m=as.data.frame(scz_meta_dge_m)
colnames(scz_meta_dge_m) = c("beta", "SE", "p")
rownames(scz_meta_dge_m) = genes
scz_meta_dge_m$fdr = p.adjust(scz_meta_dge_m$p, "fdr")
scz_meta_dge_m$symbol=SCZ_exprs_list[[1]]$annot$external_gene_id[match(genes, SCZ_exprs_list[[1]]$annot$ensembl_gene_id)]


scz_meta_dge_f=as.data.frame(scz_meta_dge_f)
colnames(scz_meta_dge_f) = c("beta", "SE", "p")
rownames(scz_meta_dge_f) = genes
scz_meta_dge_f$fdr = p.adjust(scz_meta_dge_f$p, "fdr")
scz_meta_dge_f$symbol=SCZ_exprs_list[[1]]$annot$external_gene_id[match(genes, SCZ_exprs_list[[1]]$annot$ensembl_gene_id)]

saveRDS(scz_meta_dge_m, file = "out/scz_meta_dge_m.rds")
saveRDS(scz_meta_dge_f, file = "out/scz_meta_dge_f.rds")

#correlation####
scz_meta_dge_m <- readRDS("file = out/scz_meta_dge_m.rds")
scz_meta_dge_f <- readRDS("file = out/scz_meta_dge_f.rds")
cor(scz_meta_dge_m$beta, scz_meta_dge_f$beta, use = "complete", method = "spear")
plot(scz_meta_dge_m$beta ~ scz_meta_dge_f$beta)


#BD####
files = dir("./working_data/01_unbiased/", pattern="_sex")
files = files[1:3]
n=length(files)
BD_exprs_list = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/01_unbiased/", files[[i]], sep=""))
  BD_exprs_list[[i]]$exprs = exprs 
  BD_exprs_list[[i]]$meta = meta
  BD_exprs_list[[i]]$annot = annot
  rm(exprs,meta,annot)
}

genes = intersect(rownames(BD_exprs_list[[1]]$exprs), rownames(BD_exprs_list[[2]]$exprs))
for(i in 3:n) genes = intersect(genes, rownames(BD_exprs_list[[i]]$exprs))
for(i in 1:n) BD_exprs_list[[i]]$exprs = BD_exprs_list[[i]]$exprs[match(genes,rownames(BD_exprs_list[[i]]$exprs)),]

BD_exprs_m = data.frame(row.names = genes)
BD_meta_m = data.frame(Study=NA, Subject=NA, Group=NA)
BD_exprs_f = data.frame(row.names = genes)
BD_meta_f = data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  idx <- which(BD_exprs_list[[i]]$meta[,"Sex"] == "M")
  
  BD_exprs_m = cbind(BD_exprs_m, BD_exprs_list[[i]]$exprs[,idx])
  BD_meta_m = rbind(BD_meta_m, BD_exprs_list[[i]]$meta[idx,c("Study", "Subject","Group")])
  
  BD_exprs_f = cbind(BD_exprs_f, BD_exprs_list[[i]]$exprs[,-idx])
  BD_meta_f = rbind(BD_meta_f, BD_exprs_list[[i]]$meta[-idx,c("Study", "Subject","Group")])
  
}
BD_meta_m = BD_meta_m[-1,]
to_keep = BD_meta_m$Group %in% c("CTL", "BD")
BD_exprs_m = BD_exprs_m[,to_keep]; BD_meta_m = BD_meta_m[to_keep,]
BD_meta_m$Group = factor(BD_meta_m$Group)
BD_meta_m$Group <- relevel(BD_meta_m$Group, "CTL")
BD_meta_m$Study = as.factor(BD_meta_m$Study)
BD_meta_m$Subject = as.factor(BD_meta_m$Subject)

BD_meta_f = BD_meta_f[-1,]
to_keep = BD_meta_f$Group %in% c("CTL", "BD")
BD_exprs_f = BD_exprs_f[,to_keep]; BD_meta_f = BD_meta_f[to_keep,]
BD_meta_f$Group = factor(BD_meta_f$Group)
BD_meta_f$Group <- relevel(BD_meta_f$Group, "CTL")
BD_meta_f$Study = as.factor(BD_meta_f$Study)
BD_meta_f$Subject = as.factor(BD_meta_f$Subject)

bd_meta_dge_m = matrix(NA, nrow=nrow(BD_exprs_m), ncol=3)
bd_meta_dge_f = matrix(NA, nrow=nrow(BD_exprs_f), ncol=3)
for(i in 1:nrow(BD_exprs_m)) {
  if(i%%100==0) print(i)
  expr_m = as.numeric(BD_exprs_m[i,])
  expr_f = as.numeric(BD_exprs_f[i,])
  tryCatch({
    bd_meta_dge_m[i,] = summary(lme(expr_m~ Group + Study,data = BD_meta_m, random=~1|Subject))$tTable[2,c(1,2,5)]
    bd_meta_dge_f[i,] = summary(lme(expr_f~ Group + Study,data = BD_meta_f, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e) print(e))
}

bd_meta_dge_m=as.data.frame(bd_meta_dge_m)
colnames(bd_meta_dge_m) = c("beta", "SE", "p")
rownames(bd_meta_dge_m) = genes
bd_meta_dge_m$fdr = p.adjust(bd_meta_dge_m$p, "fdr")
bd_meta_dge_m$symbol=BD_exprs_list[[1]]$annot$external_gene_id[match(genes, BD_exprs_list[[1]]$annot$ensembl_gene_id)]


bd_meta_dge_f=as.data.frame(bd_meta_dge_f)
colnames(bd_meta_dge_f) = c("beta", "SE", "p")
rownames(bd_meta_dge_f) = genes
bd_meta_dge_f$fdr = p.adjust(bd_meta_dge_f$p, "fdr")
bd_meta_dge_f$symbol=BD_exprs_list[[1]]$annot$external_gene_id[match(genes, BD_exprs_list[[1]]$annot$ensembl_gene_id)]

saveRDS(bd_meta_dge_m, file = "out/bd_meta_dge_m.rds")
saveRDS(bd_meta_dge_f, file = "out/bd_meta_dge_f.rds")

#correlation####
bd_meta_dge_m <- readRDS(file = "out/bd_meta_dge_m.rds")
bd_meta_dge_f <- readRDS(file = "out/bd_meta_dge_f.rds")
cor(bd_meta_dge_m$beta, bd_meta_dge_f$beta, use = "complete", method = "spear")
plot(bd_meta_dge_m$beta ~ bd_meta_dge_f$beta)

#BD SCZ correlation####
library(ggplot2)

bd_meta_comb <- cbind(bd_meta_dge_f, bd_meta_dge_m)
colnames(bd_meta_comb) <- c("beta_f", "SE_f", "p_f", "fdr_f", "symbol", "beta_m", "SE_m", "p_m", "fdr_m", "symbol")
bd_meta_comb$fdr_mean <- rowMeans(bd_meta_comb[,c("fdr_f", "fdr_m")], na.rm = T)

scz_meta_comb <- cbind(scz_meta_dge_f, scz_meta_dge_m)
colnames(scz_meta_comb) <- c("beta_f", "SE_f", "p_f", "fdr_f", "symbol", "beta_m", "SE_m", "p_m", "fdr_m", "symbol")
scz_meta_comb$fdr_mean <- rowMeans(scz_meta_comb[,c("fdr_f", "fdr_m")], na.rm = T)

ggplot(bd_meta_comb, aes(x = beta_f, y = beta_m, color = fdr_mean)) + geom_point(size = 1.5, alpha = .5, shape = 18) +
  geom_abline(slope = cor(bd_meta_dge_m$beta, bd_meta_dge_f$beta, use = "complete", method = "spear"), intercept = 0, color = "darkred") +
  xlab("Beta, females") + ylab("Beta, males") + ggtitle("Bipolar Disorder") +
  scale_color_gradient(name = "FDR", low = "#132B43", high = "#56B1F7") +
  scale_x_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + scale_y_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + 
  annotate("text", x = .7, y = .73, label = paste("R =", round(cor(bd_meta_dge_m$beta, bd_meta_dge_f$beta, use = "complete", method = "spear"),2)))
ggsave("img/BP_sex_cor.png")
ggplot(scz_meta_comb, aes(x = beta_f, y = beta_m, color = fdr_mean)) + geom_point(size = 1.5, alpha = .5, shape = 18) +
  geom_abline(slope = cor(scz_meta_dge_m$beta, scz_meta_dge_f$beta, use = "complete", method = "spear"), intercept = 0, color = "darkred") +
  xlab("Beta, females") + ylab("Beta, males") + ggtitle("Schizophrenia") +
  scale_color_gradient(name = "FDR", low = "#2B4313", high = "#AAB3A0") +
  scale_x_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + scale_y_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + 
  annotate("text", x = .7, y = .73, label = paste("R =", round(cor(scz_meta_dge_m$beta, scz_meta_dge_f$beta, use = "complete", method = "spear"),2)))
ggsave("img/SCZ_sex_cor.png")

cor(scz_meta_dge_m$beta, bd_meta_dge_m$beta, use = "complete", method = "spear")
plot(scz_meta_dge_m$beta ~ bd_meta_dge_m$beta)
cor(scz_meta_dge_f$beta, bd_meta_dge_f$beta, use = "complete", method = "spear")
plot(scz_meta_dge_f$beta ~ bd_meta_dge_f$beta)

scz_meta_dge <- readRDS(file = "out/scz_meta_dge.rds")
bd_meta_dge <- readRDS(file = "out/bd_meta_dge.rds")

cor(scz_meta_dge$beta, bd_meta_dge$beta, use = "complete", method = "spear")

hist(bd_meta_dge$beta)
hist(bd_meta_dge$fdr)
hist(scz_meta_dge$beta)
hist(scz_meta_dge$fdr)

scz_bd_comb <- cbind(scz_meta_dge, bd_meta_dge)
colnames(scz_bd_comb) <- c("beta_scz", "SE_scz", "p_scz", "fdr_scz", "symbol", "beta_bd", "SE_bd", "p_bd", "fdr_bd", "symbol")
scz_bd_comb$fdr_mean <- rowMeans(scz_bd_comb[,c("fdr_scz", "fdr_bd")], na.rm = T)

ggplot(scz_bd_comb, aes(x = beta_scz, y = beta_bd, color = fdr_mean)) + geom_point(size = 1.5, alpha = .5, shape = 18) +
  geom_abline(slope = cor(scz_meta_dge$beta, bd_meta_dge$beta, use = "complete", method = "spear"), intercept = 0, color = "darkred") +
  xlab("Beta, Schizophrenia") + ylab("Beta, Bipolar Disorder") + ggtitle("BD vs SCZ, all patients") +
  scale_color_gradient(name = "FDR", low = "#7E1111", high = "#CB9F9F") +
  scale_x_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + scale_y_continuous(limits = c(-.6,.8), breaks = seq(-.5,.8,.2)) + 
  annotate("text", x = .7, y = .73, label = paste("R =", round(cor(scz_meta_dge$beta, bd_meta_dge$beta, use = "complete", method = "spear"),2)))
ggsave("img/SCZ_BD_cor.png")

#save all data for permutations####
#MALE####
scz_bd_meta_dge_m = vector(mode="list", length=2)
scz_bd_meta_dge_m[[1]]$exprs = SCZ_exprs_m; scz_bd_meta_dge_m[[1]]$meta = SCZ_meta_m
scz_bd_meta_dge_m[[2]]$exprs = BD_exprs_m; scz_bd_meta_dge_m[[2]]$meta = BD_meta_m 
names(scz_bd_meta_dge_m) = c("SCZ", "BD")

genes = intersect(rownames(scz_bd_meta_dge_m[[1]]$exprs), rownames(scz_bd_meta_dge_m[[2]]$exprs))

for(i in 1:2) scz_bd_meta_dge_m[[i]]$exprs = scz_bd_meta_dge_m[[i]]$exprs[match(genes, rownames(scz_bd_meta_dge_m[[i]]$exprs)),]
saveRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge_m.RData", scz_bd_meta_dge_m)

#FEMALE####
scz_bd_meta_dge_f = vector(mode="list", length=2)
scz_bd_meta_dge_f[[1]]$exprs = SCZ_exprs_f; scz_bd_meta_dge_f[[1]]$meta = SCZ_meta_f
scz_bd_meta_dge_f[[2]]$exprs = BD_exprs_f; scz_bd_meta_dge_f[[2]]$meta = BD_meta_f 
names(scz_bd_meta_dge_f) = c("SCZ", "BD")

genes = intersect(rownames(scz_bd_meta_dge_f[[1]]$exprs), rownames(scz_bd_meta_dge_f[[2]]$exprs))

for(i in 1:2) scz_bd_meta_dge_f[[i]]$exprs = scz_bd_meta_dge_f[[i]]$exprs[match(genes, rownames(scz_bd_meta_dge_f[[i]]$exprs)),]
saveRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge_f.RData", scz_bd_meta_dge_f)

#for GO####
save(scz_bd_comb, scz_meta_comb, bd_meta_comb, 
     bd_meta_dge_m, bd_meta_dge_f, scz_meta_dge_m, scz_meta_dge_f, file = "working_data/01_unbiased/meta_for_GO.RData")
