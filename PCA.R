library(data.table)
library(DESeq2)
library(ggplot2)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Flowering/")

dds <- readRDS("counts/DESeqCasFlower.rds")

#varian-stabilizing 
vsd <- vst(dds , blind = T)

#pca on transformed data
pca <- prcomp(t(assay(vsd)), center = T, scale. = F)
##variance explained in percentages 
pca.var <- pca$sdev^2
pca.var.perc <- round(100 * pca.var / sum(pca.var), 1)

#pca dataframe 
pca.df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  time = colData(dds)$time,
  geno = colData(dds)$geno,
  condition = colData(dds)$condition,
  rep = colData(dds)$rep
)

#time color
colors <- c("T0" = "#1b9e77", "T1" = "#d95f02", "T2" = "#7570b3")

#pca plot
ggplot(pca.df, aes(x = PC1, y = PC2, color = time, shape = geno, alpha = condition)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_alpha_manual(
    values = c("Light" = 1, "NoLight" = 0.5),
    name = "Condition"   # <-- adds legend title
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", pca.var.perc[1], "%)"),
    y = paste0("PC2 (", pca.var.perc[2], "%)"),
    color = "Time point",
    shape = "Genotype"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 4, alpha = 1)),
    alpha = guide_legend(override.aes = list(size = 4))  
  )

