
library(limma)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gplots)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(biomaRt)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 30)) 
legend.title = element_text(colour = "steelblue", face = "bold.italic", family = "Helvetica")
legend.text = element_text(face = "italic", color = "steelblue4", family = "Helvetica")
axis.title = element_text(family = "Helvetica", size = (15), color = "steelblue4")
axis.text = element_text(family = "Courier", color = "cornflowerblue", size = (15))

options(stringsAsFactors = F)
setwd("E:/pro/Omics Dade/Data")
files <- list.files(".","")
cn <- lapply(files, read.delim, header = F)

cn <- do.call(cbind, cn)

rownames(cn) <- gsub("[.].*", "", cn$V1)

cn <- cn[, -seq(1, ncol(cn), 2)]

cn$gene <- rownames(cn)
genes <- cn$gene

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- cn$gene
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = genes, mart = mart)

cn <- merge(cn,G_list,by.x = "gene",by.y = "ensembl_gene_id")
cn <- cn[!duplicated(cn$hgnc_symbol),]

#plt <- distinct(cn, hgnc_symbol, .keep_all = T)
#datab <- aggregate(.~hgnc_symbol, cn, mean)

rownames(cn) <- cn$hgnc_symbol


cn <- subset(cn, select = -c(hgnc_symbol, gene))
head(cn)
dim(cn)

################################################################################################
#####################################Preparation#################################################
################################################################################################

#SN gr

gr <- c(rep("CN", 2), "PD", rep("CN", 3), rep("PD", 6), rep("CN", 2), rep("PD", 2), "CN", rep("PD", 6), "CN", "CN", "PD", "CN", "PD", "CN")
colnames(cn) <- gr
head(cn)

grT <- paste0(gr, ".", seq_along(gr))
colnames(cn) <- grT
head(cn)

##################################################################
################Differential Expression Genes#####################
##################################################################

gr <- factor(gr)


colData <- data.frame(group = gr, type = "single-end")
cds <- DESeqDataSetFromMatrix(cn, colData, design = ~group)

cds <- DESeq(cds)
cnt <- log2(1+counts(cds, normalize = T))

setwd("E:/pro/Omics Dade/Plots")


pc <- prcomp(g)
ex.scale <- t(scale(t(g), scale = F))
pc1 <- prcomp(ex.scale)
pcr <- data.frame(pc1$r[,1:3], Group = gr)
pcrt <- data.frame(pc1$r[,1:3], Group = grT)
pdf("PCgrid.pdf", width = 15)
p1 <- ggplot(pcr, aes(PC1,PC2, color =Group))+ geom_point(size = 5)+ theme_bw(base_size = 15)+ mynamestheme+ labs(title = "PCA of all samples", y = "Principal Component Analysis 2", x = "Principal Component Analysis 1")
p2 <- ggplot(pcrt, aes(PC1,PC2, label =Group))+ geom_text(size = 5)+ theme_bw(base_size = 10)+ mynamestheme+ labs(y = "Principal Component Analysis 2", x = "Principal Component Analysis 1")
p1
p2
grid.arrange(p1, p2, nrow = 2)
dev.off()

g <- subset(cnt, select = -c(CN.1, PD.16, PD.19, PD.22, PD.18, PD.12, CN.25, CN.5, PD.11, CN.24))

gr <- c(rep("CN", 1), "PD", rep("CN", 2), rep("PD", 4), rep("CN", 2), rep("PD", 1), "CN", rep("PD", 3), "PD", "CN", "PD", "CN")

grT <- colnames(g)

length(gr)
length(grT)
dim(g)

pc <- prcomp(g)
ex.scale <- t(scale(t(g), scale = F))
pc1 <- prcomp(ex.scale)
pcr <- data.frame(pc1$r[,1:3], Group = gr)
pcrt <- data.frame(pc1$r[,1:3], Group = grT)
pdf("PCgrid4.pdf", width = 15)
p1 <- ggplot(pcr, aes(PC1,PC2, color =Group))+ geom_point(size = 5)+ theme_bw(base_size = 15)+ mynamestheme+ labs(title = "PCA of all samples", y = "Principal Component Analysis 2", x = "Principal Component Analysis 1")
p2 <- ggplot(pcrt, aes(PC1,PC2, label =Group))+ geom_text(size = 5)+ theme_bw(base_size = 10)+ mynamestheme+ labs(y = "Principal Component Analysis 2", x = "Principal Component Analysis 1")
p1
p2
grid.arrange(p1, p2, nrow = 2)
dev.off()


#Remove 10 samples SN

gr <- c(rep("CN", 2), "PD", rep("CN", 3), rep("PD", 6), rep("CN", 2), rep("PD", 2), "CN", rep("PD", 6), "CN", "CN", "PD", "CN", "PD", "CN")

grT <- paste0(gr, ".", seq_along(gr))

ex <- subset(cn, select = -c(CN.1, PD.16, PD.19, PD.22, PD.18, PD.12, CN.25, CN.5, PD.11, CN.24))
gr <- c(rep("CN", 1), "PD", rep("CN", 2), rep("PD", 4), rep("CN", 2), rep("PD", 1), "CN", rep("PD", 3), "PD", "CN", "PD", "CN")

grT <- colnames(ex)
colnames(ex) <- gr
length(gr)
length(grT)
dim(ex)

##################################################################
################Differential Expression Genes#####################
##################################################################

gr <- factor(gr)


colData <- data.frame(group = gr, type = "single-end")
cds <- DESeqDataSetFromMatrix(ex, colData, design = ~group)

cds <- DESeq(cds)
cnt <- log2(1+counts(cds, normalize = T))
dif <- data.frame(results(cds, c("group", "CN", "PD"))) 

dif$padj <- p.adjust(dif$pvalue, method = "BH")
dif <- dif[order(dif$padj),]

till.up <- subset(dif, log2FoldChange > 1 & pvalue < 0.05)
till.down <- subset(dif, log2FoldChange < -1 & pvalue < 0.05)

dim(till.up)
dim(till.down)
setwd("E:/pro/Omics Dade/DEGs")
write.table(dif, "DEGs.txt", quote = F, row.names = T, sep = "\t")
write.table(cn, "Counts.txt", quote = F, row.names = T, sep = "\t")
write.table(till.up, "Up-regulated 1114 genes.txt", quote = F, row.names = T, sep = "\t")
write.table(till.down, "Down-regulated 495 genes.txt", quote = F, row.names = T, sep = "\t")


########################### Volcano Plot #############################################

library(ggrepel)
install.packages("ggrepel")

p <- ggplot(data = dif, aes(x = log2FoldChange, y = -log10(pvalue), col = log2FoldChange)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept = c(-1,1), col = "red") + geom_hline(yintercept = -log10(0.05), col = "red")
p2

#====================================================

dif$Diffexpressed <- "NO"
dif$Diffexpressed[dif$log2FoldChange > 1 & dif$pvalue < 0.05] <- "UP"
dif$Diffexpressed[dif$log2FoldChange < -1 & dif$pvalue < 0.05] <- "DOWN"

p <- ggplot(data = dif, aes(x = log2FoldChange, y = -log10(pvalue), col = Diffexpressed)) + geom_point() + theme_minimal()
p
p2 <- p + geom_vline(xintercept = c(-1,1), col = "red") + geom_hline(yintercept = -log10(0.05), col = "red")
p2

p3 <- p2 + scale_color_manual(values = c("blue", "black", "red"))
p3

setwd("E:/pro/Omics Dade/Plots")

pdf("SN_Volcano plot_sample.pdf")
p3
dev.off()





















