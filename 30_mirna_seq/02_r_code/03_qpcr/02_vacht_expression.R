#CHAT EXPRESSION IN LA-N-2 CNTF DIFFERENTIATION####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(ggplot2)

expr <- readRDS("./raw_data/raw_expression_qpcr_la2_vacht.rds")

#pval
t.test(expr$Expression[expr$Group == "II_CNTF"], 
       expr$Expression[expr$Group == "II_CON"])
t.test(expr$Expression[expr$Group == "IV_CNTF"], 
       expr$Expression[expr$Group == "IV_CON"])
