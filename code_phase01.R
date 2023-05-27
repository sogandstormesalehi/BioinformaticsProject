setwd("C:\\Users\\vcc\\Desktop\\bioproj")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rtsne)
library(MASS)
library(dplyr)
library(M3C)

#download data
series <- "GSE48558"
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "Data\\")
gset <- gset[[1]]

#set groupings, X is for data we don't need to process
gr <-  c('AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Healthy', 'X', 'X', 'X', 'Healthy', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Healthy', 'X', 'Healthy', 'X', 'X', 'X', 'Healthy', 'X', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'X', 'Healthy', 'X', 'X', 'Healthy', 'Healthy', 'X', 'X', 'Healthy', 'Healthy', 'X', 'Healthy', 'X', 'Healthy', 'X', 'Healthy', 'X', 'Healthy', 'X', 'Healthy', 'X', 'X', 'X', 'Healthy', 'X', 'X', 'X', 'Healthy', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'AML', 'AML', 'Healthy', 'AML', 'AML', 'AML', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')

#get expression matrix from gset
ex <- exprs(gset)

#draw boxplot to see if the data is normalized and log scaled
pdf("Results\\boxplot.pdf", width=170, height=44)
boxplot(ex, col="mistyrose")
dev.off()

#draw heatmap to get a sense of the correlation between different groups
pdf("Results\\CorHeatmap.pdf", width=60, height=60)
pheatmap(cor(ex), labels_row = gr, labels_col = gr, color=cm.colors(256), border_color = NA)
dev.off()

#mean center the matrix for PC analysis
ex.scale <- t(scale(t(ex), scale=F))

#dimensionality reduction using PCA
pc <- prcomp(ex.scale)

#draw plots to see the variace along different PCs
pdf("Results\\PCA.pdf")
plot(pc, col="#CC99FF")
plot(pc$x[,1:2], col="#6699FF")
dev.off()

#show samples in the new, reduced space
pcr <- data.frame(pc$rotation[,1:3], Group=gr)
pdf("Results\\PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 2) + theme_bw()
dev.off()

#dimensionality reduction using t_SNE
ex.frame <- data.frame(ex.scale)
TSNE <- tsne(ex.frame, perplex=44)

#show samples in the new, reduced space
tsne.frame <- data.frame(TSNE$data, Group = gr)
pdf("Results\\tSNE.pdf")
ggplot(tsne.frame, aes(X1, X2, color = Group)) + geom_point(size = 2) + theme_bw()
dev.off()

#dimensionality reduction using MDS
ex.frame2 <- data.frame(t(ex.scale))
mds <- ex.frame2  %>%
  dist('euclidean') %>%          
  isoMDS() %>%
  .$points %>%
  as_tibble()

#show samples in the new, reduced space
colnames(mds) <- c("X", "Y")
mds.frame <- data.frame(X = mds$X, Y = mds$Y, Group = gr)
pdf("Results\\MDS.pdf")
ggplot(mds.frame, aes(X, Y, color = Group)) + geom_point(size = 2) + theme_bw()
dev.off()


#new groupings based on the normal samples' source name
gr.second <- c('AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Granulocytes', 'X', 'X', 'X', 'Granulocytes', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'B cells', 'X', 'T cells', 'X', 'X', 'X', 'Granulocytes', 'X', 'Granulocytes', 'Monocytes', 'Monocytes', 'B cells', 'X', 'T cells', 'X', 'X', 'T cells', 'T cells', 'X', 'X', 'T cells', 'T cells', 'X', 'B cells', 'X', 'T cells', 'X', 'B cells', 'X', 'T cells', 'X', 'CD34+HSPC', 'X', 'X', 'X', 'CD34+HSPC', 'X', 'X', 'X', 'CD34+HSPC', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'AML', 'AML', 'T cells', 'AML', 'AML', 'AML', 'B cells', 'B cells', 'B cells', 'B cells', 'B cells', 'B cells', 'B cells', 'T cells', 'Monocytes', 'Monocytes', 'Monocytes', 'Monocytes', 'Granulocytes', 'T cells', 'T cells', 'T cells', 'T cells', 'T cells', 'T cells', 'T cells')

#filter out the X samples
ex.filtered <- ex[,gr.second != 'X']
gr.filtered <- gr.second[gr.second != 'X']

#draw heatmap to compare different groups' correlation
pdf("Results\\CorHeatmap_second.pdf", width=80, height=80)
pheatmap(cor(ex.filtered), labels_row = gr.filtered, labels_col = gr.filtered, color=cm.colors(256), border_color = NA)
dev.off()

#draw tsne plot considering the new grouping
tsne.frame2 <- data.frame(TSNE$data, Group = gr.second)
pdf("Results\\tSNE_second.pdf")
ggplot(tsne.frame2, aes(X1, X2, color = Group)) + geom_point(size = 2) + theme_bw()
dev.off()
