#SUPPLEMENT TO MULTI-SHEET XLSX DataS####
home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(rJava)
library(xlsx)

#File S1####
#meta GO
file.remove("out/DataS1.csv")
#master
M <- read.table("out/sheet1.csv", header = T)
M$set <- "sex-independent"
# write.xlsx(M, file = "out/DataS1.xlsx", sheetName = "SCZvsBD", append = T, row.names = F)

#disease differences per sex A + B
#A
#left
A1 <- read.table("out/sheet5.csv", header = T)
A1$set <- "SCZ-biased F"
# write.xlsx(A1, file = "out/DataS1.xlsx", sheetName = "SCZ-biased F", append = T, row.names = F)
#right
A2 <- read.table("out/sheet4.csv", header = T)
A2$set <- "SCZ-biased M"
# write.xlsx(A2, file = "out/DataS1.xlsx", sheetName = "SCZ-biased M", append = T, row.names = F)
#B
#left
B1 <- read.table("out/sheet7.csv", header = T)
B1$set <- "BD-biased F"
# write.xlsx(B1, file = "out/DataS1.xlsx", sheetName = "BD-biased F", append = T, row.names = F)
#right
B2 <- read.table("out/sheet6.csv", header = T)
B2$set <- "BD-biased M"
# write.xlsx(B2, file = "out/DataS1.xlsx", sheetName = "BD-biased M", append = T, row.names = F)

#sex differences per disease C + D
#C
#left
C1 <- read.table("out/sheet10.csv", header = T)
C1$set <- "M-biased BD"
# write.xlsx(C1, file = "out/DataS1.xlsx", sheetName = "M-biased BD", append = T, row.names = F)
#right
C2 <- read.table("out/sheet13.csv", header = T)
C2$set <- "M-biased SCZ"
# write.xlsx(C2, file = "out/DataS1.xlsx", sheetName = "M-biased SCZ", append = T, row.names = F)
#D
#left
D1 <- read.table("out/sheet9.csv", header = T)
D1$set <- "F-biased BD"
# write.xlsx(D1, file = "out/DataS1.xlsx", sheetName = "F-biased BD", append = T, row.names = F)
#right
D2 <- read.table("out/sheet12.csv", header = T)
D2$set <- "F-biased SCZ"
# write.xlsx(D2, file = "out/DataS1.xlsx", sheetName = "F-biased SCZ", append = T, row.names = F)

S1 <- rbind(M[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            A1[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            A2[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            B1[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            B2[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            C1[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            C2[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            D1[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")], 
            D2[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "weight", "set")])

colnames(S1) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "adj.p.val", "set")
write.table(S1, file = "out/DataS1.csv", sep = ";", row.names = F, quote = T)

#File S2####
file.remove("out/DataS2.csv")
#permutation mirs, tfs, tf mirs (Fig3)
primate_specific_mirs <- readRDS(file = "./working_data/single_cell/primate_specific_mirs.rds")
non_ps_mirs <- readRDS(file = "./working_data/single_cell/non_primate_specific_mirs.rds")
tfs <- readRDS(file = "./working_data/single_cell/significant_tfs.rds")
tf_mirs <- readRDS(file = "./working_data/single_cell/significant_tf_mirs.rds")

colnames(primate_specific_mirs) <- c("name", "number_of_targets", "perm.p.val")
colnames(non_ps_mirs) <- c("name", "number_of_targets", "perm.p.val")
colnames(tfs) <- c("name", "number_of_targets", "perm.p.val")
colnames(tf_mirs) <- c("name", "number_of_targets", "perm.p.val")

primate_specific_mirs$biotype <- "primate-specific_miR"
non_ps_mirs$biotype <- "conserved_miR"
tfs$biotype <- "transcription_factor"
tf_mirs$biotype <- "TF-targeting_miR"

perm <- rbind(primate_specific_mirs,
              non_ps_mirs,
              tfs,
              tf_mirs)

write.table(perm, file = "out/DataS2.csv", row.names = F, quote = F, sep = ";")

#File S3####
file.remove("out/DataS3.csv")
#de mirs la2/5 (Fig4)
lnames <- load("working_data/mirna_seq/LA2_LA5_DESeq_results_new_seq.RData")
lnames <- lnames[-grep("vs", lnames)] #no interaction terms
de_l <- list()
for(nam in lnames){
  de_l[[nam]] <- data.frame(get(nam))
}
de_f <- plyr::ldply(de_l)
colnames(de_f)[1] <- "condition"

write.table(de_f, file = "out/DataS3.csv", row.names = F, quote = F, sep = ";")

#File S4####
#literature mirs (Fig4), created in 30/02/02/03_

#File S5####
file.remove("out/DataS5.csv")
#miR family GO
mir_fam_go_list <- readRDS(file = "working_data/mir_family_GO_enrich_list.rds")
go <- plyr::ldply(mir_fam_go_list)
go <- go[, c(".id","GO.ID","Term","Annotated","Significant","Expected","weight")]
colnames(go) <- c("Family_name","GO.ID","Term","Annotated","Significant","Expected","adj.p.val")

write.table(go, file = "out/DataS5.csv", row.names = F, quote = T, sep = ";")

#File S6####
#network, created in 30/02/04/02_

#File S7####
#filtered genes (Fig6), created in 60/02a_
#File S8####
#filtered mirs (Fig6), created in 30/02/06/04_