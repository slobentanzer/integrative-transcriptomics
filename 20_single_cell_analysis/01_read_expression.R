#COMPILE SINGLE CELL SEQ DATA####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(RNeo4j)
graph <- startGraph("http://localhost:7474/db/data/")

darmanis <- readRDS("./raw_data/darmanis_chol_exprs.rds")
habib <- readRDS("./raw_data/habib_chol_exprs.rds")
tasic <- readRDS("./raw_data/tasic_chol_exprs.rds")
zeisel <- readRDS("./raw_data/zeisel_chol_exprs.rds")

#normalize
darmanis[,1:9] <- apply(darmanis[,1:9], 2, FUN = function(X) (X - min(X))/diff(range(X)))
habib[,1:9] <- apply(habib[,1:9], 2, FUN = function(X) (X - min(X))/diff(range(X)))
tasic[,1:4] <- apply(tasic[,1:4], 2, FUN = function(X) (X - min(X))/diff(range(X)))
zeisel[,1:4] <- apply(zeisel[,1:4], 2, FUN = function(X) (X - min(X))/diff(range(X)))

#heatmap
library(svglite)
library(pheatmap)
svg("./img/darmanis_heatmap.svg", height = 10)
pheatmap(as.matrix(darmanis[,1:9]),
         cellwidth = 10,
         cellheight = 10,
         border_color = NA#,
         # filename = "darmanis_heatmap.jpg"
         )
dev.off()
svg("./img/habib_heatmap.svg", height = 10)
pheatmap(as.matrix(habib[,1:9]),
         cellwidth = 10, 
         cellheight = 10,
         border_color = NA#,
         # filename = "habib_heatmap.jpg"
)
dev.off()
svg("./img/tasic_heatmap.svg", height = 10)
pheatmap(as.matrix(tasic[,1:4]),
         cellwidth = 10, 
         cellheight = 10,
         border_color = NA#,
         # filename = "tasic_heatmap.jpg"
)
dev.off()
svg("./img/zeisel_heatmap.svg", height = 10)
pheatmap(as.matrix(zeisel[,1:4]),
         cellwidth = 10, 
         cellheight = 10,
         border_color = NA#,
         # filename = "zeisel_heatmap.jpg"
         )
dev.off()

#boxplot
library(reshape2)
darmanis.melt <- melt(darmanis[,!colnames(darmanis) %in% c("group", "order")])
habib.melt <- melt(habib)
tasic.melt <- melt(tasic)
zeisel.melt <- melt(zeisel)

comp <- rbind(darmanis.melt, habib.melt, tasic.melt, zeisel.melt)
library(RColorBrewer)
col <- brewer.pal(n = 8, name = "Paired")
col1 <- col[c(1,3,5,7)]
col2 <- col[c(2,4,6,8)]
ggplot(comp, aes(name, value, fill = Study, color = Study)) + 
  ylab("Log-fold expression") + xlab("Gene symbol") +
  theme_minimal() + #coord_flip() +
  scale_color_manual(values = col2) + 
  scale_fill_manual(values = col1) + 
  #geom_hline(yintercept = 0) +
  geom_boxplot(outlier.size = .3, lwd = .2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("single_cell_comp_boxplot.png", width = 7, height = 4, dpi = 1000)
