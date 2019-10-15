#DE MIR FAMILIES####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

options(stringsAsFactors = F)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

#mirbase data####
library(readr)
miFam <- read.table("raw_data/miFam.dat", fill = T, sep = " ")
miFam <- data.frame(miFam)
ac_idx <- grep("AC", miFam$V1)
id_idx <- grep("ID", miFam$V1)
hsa_idx <- grep("hsa", miFam$V6)
idx <- sort(c(ac_idx, id_idx, hsa_idx))
miFam <- miFam[idx,]

#de data####
mir.matrix <- readRDS(file = "working_data/mir_de_matrix_countchange.rds")

de_mirs <- rownames(mir.matrix)[rowSums(mir.matrix) > 0]

#find de in mirbase####
de_mirs <- gsub("-[3-5]p", "", de_mirs)
de_mirs <- tolower(de_mirs)
mean(de_mirs %in% miFam$V6)
de_mirs[!de_mirs %in% miFam$V6]

#..-1/-2 etc inconsistent. remove####
miFam$V7 <- miFam$V6
miFam$V7[grep("-[1-9]$", miFam$V7)] <- unlist(lapply(strsplit(miFam$V6[grep("-[1-9]$", miFam$V6)], "-"), function(x) paste(x[1:3], collapse = "-")))
de_mirs[grep("-[1-9]$", de_mirs)] <- unlist(lapply(strsplit(de_mirs[grep("-[1-9]$", de_mirs)], "-"), function(x) paste(x[1:3], collapse = "-")))
mean(de_mirs %in% miFam$V6)
mean(de_mirs %in% miFam$V7)
de_mirs[!de_mirs %in% miFam$V7] #17 left, not high confidence

#miFam to list####
head(miFam)
miFam_list <- vector(mode = "list", length = length(grep("AC", miFam$V1)))
names(miFam_list) <- miFam$V4[grep("AC", miFam$V1)+1]
for(i in 1:nrow(miFam)){
  if(miFam$V1[i] == "AC") {
    #new entry
    ac <- miFam$V4[i]
    nam <- miFam$V4[i+1]
    temp <- data.frame()
  }
  if(miFam$V1[i] == "MI") {
    #add miR
    temp <- rbind(temp, data.frame(accession = miFam$V4[i], name = miFam$V6[i], short_name = miFam$V7[i], family_accession = ac, family_name = nam))
  }
  miFam_list[[nam]] <- temp
}
miFam_list <- miFam_list[unlist(lapply(miFam_list, function(x) nrow(x)>1))]
miFam_list[["mir-10"]]

#multiple accessions for single miRs
miFam_list.red <- lapply(miFam_list, function(x) unique(x[, c("name", "short_name", "family_name", "family_accession")]))
miFam_list.red <- miFam_list.red[unlist(lapply(miFam_list.red, function(x) nrow(x)>1))]

miFam <- plyr::ldply(miFam_list.red)
miFam$.id <- NULL
miFam <- miFam[order(miFam$family_name, miFam$name),]
table(miFam$family_name)

#add mature miRs####
anc <- data.frame()
for(nam in unique(miFam$short_name)) {
  anc <- rbind(anc,
                cypher(graph, paste0("MATCH (a:ANCESTOR)--(m:MIR {species: 'HSA'})
              WHERE a.name =~ '", nam, ".*' ", 
                       "RETURN a.name, m.name"), nam = nam))
  
}

mean(is.na(match(miFam$name, anc$a.name)))
match(miFam$name, anc$a.name)
match(anc$a.name, miFam$name)

miFam <- cbind(anc, miFam[match(anc$a.name, miFam$name),])
miFam <- miFam[!is.na(miFam$name),]
miFam <- miFam[order(miFam$family_name, miFam$name, miFam$m.name),]

saveRDS(miFam, "working_data/mir_families_mirbase_v21.rds")
dplyr::count(miFam, family_name, sort = T)
miFam_list[["mir-548"]]

#DE mir family distribution####
mean(rownames(mir.matrix) %in% miFam$m.name) #about 50% are in families
rownames(mir.matrix)[!rownames(mir.matrix) %in% miFam$m.name]
