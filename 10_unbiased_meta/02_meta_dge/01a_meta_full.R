#FULL META-ANALYSIS OF BD AND SCZ SAMPLES, REPLICATING GANDAL ET AL####

rm(list=ls())
options(stringsAsFactors=F)

library(nlme)

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

#SCZ####
files = dir("./working_data/01_unbiased/", pattern="_normalized")
files = files[-grep("sex",files)]
n=length(files)
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

SCZ_exprs = data.frame(row.names = genes)
SCZ_meta = data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  SCZ_exprs = cbind(SCZ_exprs, SCZ_exprs_list[[i]]$exprs)
  SCZ_meta = rbind(SCZ_meta, SCZ_exprs_list[[i]]$meta[,c("Study", "Subject","Group")])
}
SCZ_meta = SCZ_meta[-1,]
to_keep = SCZ_meta$Group %in% c("CTL", "SCZ")
SCZ_exprs = SCZ_exprs[,to_keep]; SCZ_meta = SCZ_meta[to_keep,]
SCZ_meta$Group = factor(SCZ_meta$Group)
SCZ_meta$Study = as.factor(SCZ_meta$Study)
SCZ_meta$Subject = as.factor(SCZ_meta$Subject)
nrow(SCZ_meta[SCZ_meta$Group == "SCZ",]) #139 patients

scz_meta_dge = matrix(NA, nrow=nrow(SCZ_exprs), ncol=3)
for(i in 1:nrow(SCZ_exprs)) {
  if(i%%100==0) print(i)
  expr = as.numeric(SCZ_exprs[i,])
  tryCatch({
    scz_meta_dge[i,] = summary(lme(expr~ Group + Study,data = SCZ_meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){print(e)})
}
  
scz_meta_dge=as.data.frame(scz_meta_dge)
colnames(scz_meta_dge) = c("beta", "SE", "p")
rownames(scz_meta_dge) = genes
scz_meta_dge$fdr = p.adjust(scz_meta_dge$p, "fdr")
scz_meta_dge$symbol=SCZ_exprs_list[[1]]$annot$external_gene_id[match(genes, SCZ_exprs_list[[1]]$annot$ensembl_gene_id)]

saveRDS(scz_meta_dge, file = "out/scz_meta_dge.rds")

#BD####
files = dir("./working_data/01_unbiased/", pattern="_normalized")
files = files[-grep("sex",files)]
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

BD_exprs = data.frame(row.names = genes)
BD_meta = data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  BD_exprs = cbind(BD_exprs, BD_exprs_list[[i]]$exprs)
  #if("Group.BD" %in% colnames(BD_exprs_list[[i]]$meta)) BD_exprs_list[[i]]$meta$Group = BD_exprs_list[[i]]$meta$Group.BD ## --> only for unique controls
  BD_meta = rbind(BD_meta, BD_exprs_list[[i]]$meta[,c("Study", "Subject","Group")])
}
BD_meta = BD_meta[-1,]
to_keep = BD_meta$Group %in% c("CTL", "BD")
BD_exprs = BD_exprs[,to_keep]; BD_meta = BD_meta[to_keep,]
BD_meta$Group = factor(BD_meta$Group, levels=c("CTL", "BD"))
BD_meta$Study = as.factor(BD_meta$Study)
BD_meta$Subject = as.factor(BD_meta$Subject)
nrow(BD_meta[BD_meta$Group == "BD",]) #82 patients

bd_meta_dge = matrix(NA, nrow=nrow(BD_exprs), ncol=3)
for(i in 1:nrow(BD_exprs)) {
  if(i%%100==0) print(i)
  expr = as.numeric(BD_exprs[i,])
  tryCatch({
    bd_meta_dge[i,] = summary(lme(expr~ Group + Study,data = BD_meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta_dge=as.data.frame(bd_meta_dge)
colnames(bd_meta_dge) = c("beta", "SE", "p")
rownames(bd_meta_dge) = genes
bd_meta_dge$fdr = p.adjust(bd_meta_dge$p, "fdr")
bd_meta_dge$symbol=BD_exprs_list[[1]]$annot$external_gene_id[match(genes, BD_exprs_list[[1]]$annot$ensembl_gene_id)]

saveRDS(bd_meta_dge, file = "out/bd_meta_dge.rds")

#correlation####
bd_meta_dge <- readRDS(file = "out/bd_meta_dge.rds")
scz_meta_dge <- readRDS(file = "out/scz_meta_dge.rds")
cor(scz_meta_dge$beta, bd_meta_dge$beta, use = "complete", method = "spear")
plot(scz_meta_dge$beta ~ bd_meta_dge$beta)

#significant
nrow(bd_meta_dge)
bd_meta_dge$fdr[is.na(bd_meta_dge$fdr)] <- 1
mean(bd_meta_dge$fdr < .05) #4.2%
sum(bd_meta_dge$fdr < .05) #515

nrow(scz_meta_dge)
scz_meta_dge$fdr[is.na(scz_meta_dge$fdr)] <- 1
mean(scz_meta_dge$fdr < .05) #16.4%
sum(scz_meta_dge$fdr < .05) #2036

#both
sum(bd_meta_dge$symbol[bd_meta_dge$fdr < .05] %in% scz_meta_dge$symbol[scz_meta_dge$fdr < .05]) #376
sum(bd_meta_dge$symbol[bd_meta_dge$fdr < .05] %in% scz_meta_dge$symbol[scz_meta_dge$fdr < .05])/nrow(bd_meta_dge)*100 #3.0%

#save all data for permutations
scz_bd_meta_dge = vector(mode="list", length=2)
scz_bd_meta_dge[[1]]$exprs = SCZ_exprs; scz_bd_meta_dge[[1]]$meta = SCZ_meta 
scz_bd_meta_dge[[2]]$exprs = BD_exprs; scz_bd_meta_dge[[2]]$meta = BD_meta 
names(scz_bd_meta_dge) = c("SCZ", "BD")

genes = intersect(rownames(scz_bd_meta_dge[[1]]$exprs), rownames(scz_bd_meta_dge[[2]]$exprs))

for(i in 1:2) scz_bd_meta_dge[[i]]$exprs = scz_bd_meta_dge[[i]]$exprs[match(genes, rownames(scz_bd_meta_dge[[i]]$exprs)),]
saveRDS(file = "./working_data/01_unbiased/scz_bd_meta_dge.rds", scz_bd_meta_dge)
