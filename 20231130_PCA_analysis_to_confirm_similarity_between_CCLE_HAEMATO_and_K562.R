# this script is to investigate whether "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" is similar with K562 by PCA analysis
# 2023/11/30 made

# make directory
setwd("C:/Rdata")
dir.create("20231130_PCA_analysis_to_confirm_similarity_between_CCLE_HAEMATO_and_K562")

# activate package
library(ggplot2)
library(tibble)
library(umap)

# import table of CCLE miRNA expression
# this table is located at "\\fsz-p21.naist.jp\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
miRNA.exp.table <-read.table("CCLE_miRNA_20181103.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# extract name of cell line
cell.line <-colnames(miRNA.exp.table[,-1])

# transform table and arrange it
miRNA.exp.table <-as.data.frame(t(miRNA.exp.table))
colnames(miRNA.exp.table) <-miRNA.exp.table[1,]
rownames(miRNA.exp.table) <-NULL
miRNA.exp.table <-miRNA.exp.table[-1,]

# convert to numeric (due to transform)
miRNA.exp.table <-apply(miRNA.exp.table, 2, as.numeric)

# log2 transform
miRNA.exp.table <-as.data.frame(apply(miRNA.exp.table, 2, log2))

# name rows
rownames(miRNA.exp.table) <-cell.line

# perform PCA analysis
rpca=prcomp(miRNA.exp.table,scale=F)
pca.out <-as.data.frame(rpca$x)

# calculate contribution rate
v <-rpca$sdev^2
v <-v/sum(v)*100

# annotate "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "OTHER" and "K562"
hemato <-grep("HAEMATO",rownames(pca.out))
K562 <-grep("K562",rownames(pca.out))
pca.out[,735] <-1
pca.out[hemato,735] <-"HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
pca.out[-hemato,735] <-"OTHER"
pca.out[K562,735] <-"K562"
colnames(pca.out)[735] <-"site_primary"

# draw PCA plot
setwd("C:/Rdata/20231130_PCA_analysis_to_confirm_similarity_between_CCLE_HAEMATO_and_K562")
pca.miRNA <- ggplot(pca.out, aes(x = PC1, y = PC2, color = site_primary))
pca.miRNA  <- pca.miRNA + geom_point()+labs(x=paste0("PC1 (",signif(v[1],3),"%)"),y=paste0("PC2 (",signif(v[2],3),"%)"))
plot(pca.miRNA)
ggsave(filename = "PCA_plot_about_CCLE_miRNA_expression.pdf",plot=pca.miRNA)

# perform umap
miRNA_data_umap <- umap(miRNA.exp.table, n_components = 2)

# make table to draw umap plot 
miRNA_data_gg <- tibble(species = pca.out$site_primary, 
                             pc1 = miRNA_data_umap$layout[,1],
                             pc2 = miRNA_data_umap$layout[,2])

# draw umap plot
umap_miRNA <-ggplot(data = miRNA_data_gg, aes(x = pc1, y = pc2)) + 
  geom_point(aes(col = species), alpha = 0.6, size = 2) + 
  theme_bw()+
  theme(legend.position = "top", 
        legend.background = element_rect(colour = "gray"))

ggsave(filename = "umap_plot_about_CCLE_miRNA_expression.pdf",plot=umap_miRNA)

# import table CCLE gene expression
# this table is located at "\\fsz-p21.naist.jp\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
gene.exp.table <-read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)
gene.exp.table <-gene.exp.table[,-2]

# extract name of cell line
cell.line2 <-colnames(gene.exp.table[,-1])

# transform table and arrange it
gene.exp.table <-as.data.frame(t(gene.exp.table))
colnames(gene.exp.table) <-gene.exp.table[1,]
rownames(gene.exp.table) <-NULL
gene.exp.table <-gene.exp.table[-1,]

# convert to numeric (due to transform)
gene.exp.table <-apply(gene.exp.table, 2, as.numeric)

# add 1 TPM each column (for zero expression)
for (i in 1:ncol(gene.exp.table)) {
  gene.exp.table[,i] <-gene.exp.table[,i]+1
}

# log2 transform
gene.exp.table <-as.data.frame(apply(gene.exp.table, 2, log2))

# name rows
rownames(gene.exp.table) <-cell.line2

# perform PCA analysis
rpca=prcomp(gene.exp.table,scale=F)
pca.out <-as.data.frame(rpca$x)

# calculate contribution rate
v <-rpca$sdev^2
v <-v/sum(v)*100

# annotate "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "OTHER" and "K562"
hemato <-grep("HAEMATO",rownames(pca.out))
K562 <-grep("K562",rownames(pca.out))
pca.out[,1020] <-1
pca.out[hemato,1020] <-"HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
pca.out[-hemato,1020] <-"OTHER"
pca.out[K562,1020] <-"K562"
colnames(pca.out)[1020] <-"site_primary"

# draw PCA plot
setwd("C:/Rdata/20231130_PCA_analysis_to_confirm_similarity_between_CCLE_HAEMATO_and_K562")
pca.gene <- ggplot(pca.out, aes(x = PC1, y = PC2, color = site_primary))
pca.gene  <- pca.gene  + geom_point()+labs(x=paste0("PC1 (",signif(v[1],3),"%)"),y=paste0("PC2 (",signif(v[2],3),"%)"))
plot(pca.gene )
ggsave(filename = "PCA_plot_about_CCLE_gene_expression.pdf",plot=pca.gene)

# perform umap
gene_data_umap <- umap(gene.exp.table, n_components = 2)

# make table to draw umap plot 
gene_data_gg <- tibble(species = pca.out$site_primary, 
                             pc1 = gene_data_umap$layout[,1],
                             pc2 = gene_data_umap$layout[,2])

# draw umap plot
umap_gene <-ggplot(data = gene_data_gg, aes(x = pc1, y = pc2)) + 
  geom_point(aes(col = species), alpha = 0.6, size = 2) + 
  theme_bw()+
  theme(legend.position = "top", 
        legend.background = element_rect(colour = "gray"))

plot(umap_gene)
ggsave(filename = "umap_plot_about_CCLE_gene_expression.pdf",plot=umap_gene)
