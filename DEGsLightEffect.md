

```{r setup, include=F}

library(data.table)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(VennDiagram)
library(grid)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Flowering/")

#load - design - ~ condition * geno * time
dds <- readRDS("counts/DESeqCasFlower.rds")


#save paramters
alpha <- 0.05
lfc   <- 1
timepoints <- c("T0","T1","T2")
```


```{r}
## available coefficients
resultsNames(dds)
```


## DSC272_early - Light vs NoLight - timepoints (T0,T1,T2)

```{r}
res.DSC272.T0 <- results(dds, name = "condition_Light_vs_NoLight")

res.DSC272.T1 <- results(dds, contrast = list(c(
  "condition_Light_vs_NoLight",
  "conditionLight.timeT1"
)))

res.DSC272.T2 <- results(dds, contrast = list(c(
  "condition_Light_vs_NoLight",
  "conditionLight.timeT2"
)))
```


```{r}

#combine time points results 

res.DSC272.all <- data.table(
  gene = rownames(res.DSC272.T0),

  log2FC_T0 = res.DSC272.T0$log2FoldChange,
  padj_T0   = res.DSC272.T0$padj,

  log2FC_T1 = res.DSC272.T1$log2FoldChange,
  padj_T1   = res.DSC272.T1$padj,

  log2FC_T2 = res.DSC272.T2$log2FoldChange,
  padj_T2   = res.DSC272.T2$padj
)

head(res.DSC272.all)

```

```{r}
#counts 
# 
count.T0 <- sum(res.DSC272.all$padj_T0 < alpha & abs(res.DSC272.all$log2FC_T0) >= lfc, na.rm = TRUE)
count.T1 <- sum(res.DSC272.all$padj_T1 < alpha & abs(res.DSC272.all$log2FC_T1) >= lfc, na.rm = TRUE)
count.T2 <- sum(res.DSC272.all$padj_T2 < alpha & abs(res.DSC272.all$log2FC_T2) >= lfc, na.rm = TRUE)

# 
deg.counts.DSC272 <- data.table(
  Time = c("T0", "T1", "T2"),
  N_DEGs = c(count.T0, count.T1, count.T2)
)

deg.counts.DSC272
```

```{r}
sig.genes.DSC272 <- res.DSC272.all[(padj_T0 < alpha & abs(log2FC_T0) >= lfc) |
                                     (padj_T1 < alpha & abs(log2FC_T1) >= lfc) |
                                     (padj_T2 < alpha & abs(log2FC_T2) >= lfc)]

gene.list.DSC272 <- sig.genes.DSC272$gene


```


### Heatmap of significant genes 
```{r}
#subset for heatmap plotting 
dds.DSC272 <- dds[, colData(dds)$geno == "DSC272_early"]

#exp matrix for heatmap
vsd.DSC272 <- vst(dds.DSC272, blind = T)
mat.DSC272 <- assay(vsd.DSC272)

##keep significant genes
keep.DSC272 <- intersect(gene.list.DSC272, rownames(mat.DSC272))
mat.DSC272  <- mat.DSC272[keep.DSC272, , drop=FALSE]

#try scaled - 
#mat.scaled.DSC272 <- t(scale(t(mat.DSC272))) - no difference on the heatmap

#annotation columnes 
anno.col.DSC272 <- data.frame(
  Time      = colData(dds.DSC272)$time,
  Condition = colData(dds.DSC272)$condition
)
rownames(anno.col.DSC272) <- colnames(dds.DSC272)

#order
sample.order.DSC272 <- order(colData(dds.DSC272)$time,
                             colData(dds.DSC272)$condition)



mat.DSC272 <- mat.DSC272[, sample.order.DSC272, drop=FALSE]

anno.col.DSC272 <- anno.col.DSC272[sample.order.DSC272, , drop=FALSE]

## heatmap (all significant genes)
pheatmap(mat.DSC272,
         scale = "row",
         annotation_col = anno.col.DSC272,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "DSC272_early: sig genes")

```


### Clustering genes based on expression pattern - Fig 1.b

```{r}
#mat.DSC272

#scale 
mat.scaled.DSC272 <- t(scale(t(mat.DSC272)))

#check if there are NAs
any(is.na(mat.scaled.DSC272))
any(is.na(mat.DSC272))

##coldata
coldata.DSC272 <- as.data.table(colData(dds.DSC272), keep.rownames = "sample")

sample.order.DSC272 <- order(coldata.DSC272$time, coldata.DSC272$condition)

coldata.DSC272 <- coldata.DSC272[sample.order.DSC272]

#k-means clustering 
set.seed(123)
km.DSC272 <- kmeans(mat.scaled.DSC272, centers = 6)

#stable clustering
centers.DSC272 <- km.DSC272$centers
#order by mean expression 
cluster.order.DSC272 <- order(rowMeans(centers.DSC272))
#map
new.cluster.labels.DSC272 <- match(km.DSC272$cluster, cluster.order.DSC272)

#cluster dataframe
cluster.df.DSC272 <- data.table(
  gene = rownames(mat.scaled.DSC272),
  cluster = factor(new.cluster.labels.DSC272, levels = 1:6)
)

#order matrix by cluster
ord.DSC272 <- order(cluster.df.DSC272$cluster)
mat.scaled.DSC272 <- mat.scaled.DSC272[ord.DSC272, ]
cluster.df.DSC272 <- cluster.df.DSC272[ord.DSC272]

#row annotata
anno.row.DSC272 <- data.frame(Cluster = cluster.df.DSC272$cluster)
rownames(anno.row.DSC272) <- cluster.df.DSC272$gene


#heatmap
pheatmap(mat.scaled.DSC272,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anno.col.DSC272,
         annotation_row =anno.row.DSC272,
         show_rownames = F,
         main = "DSC272",
         cutree_rows = 6)

```



## DSC120_late  - Light vs NoLight - timepoints (T0,T1,T2)

```{r}
res.DSC120.T0 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC120_late"
                         )))

res.DSC120.T1 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC120_late",
                           "conditionLight.timeT1",
                           "conditionLight.genoDSC120_late.timeT1"
                         )))

res.DSC120.T2 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC120_late",
                           "conditionLight.timeT2",
                           "conditionLight.genoDSC120_late.timeT2"
                         )))
```





```{r}
##combine
res.DSC120.all <- data.table(
  gene = rownames(res.DSC120.T0),
  
  log2FC_T0   = res.DSC120.T0$log2FoldChange,
  padj_T0     = res.DSC120.T0$padj,
  
  log2FC_T1   = res.DSC120.T1$log2FoldChange,
  padj_T1     = res.DSC120.T1$padj,
  
  log2FC_T2   = res.DSC120.T2$log2FoldChange,
  padj_T2     = res.DSC120.T2$padj
)

```


```{r}
#counts 
count.T0 <- sum(res.DSC120.all$padj_T0 < alpha & abs(res.DSC120.all$log2FC_T0) >= lfc, na.rm = TRUE)
count.T1 <- sum(res.DSC120.all$padj_T1 < alpha & abs(res.DSC120.all$log2FC_T1) >= lfc, na.rm = TRUE)
count.T2 <- sum(res.DSC120.all$padj_T2 < alpha & abs(res.DSC120.all$log2FC_T2) >= lfc, na.rm = TRUE)

# 
deg.counts.DSC120 <- data.table(
  Time = c("T0", "T1", "T2"),
  N_DEGs = c(count.T0, count.T1, count.T2)
)

deg.counts.DSC120
```


```{r}
sig.genes.DSC120 <- res.DSC120.all[
  (padj_T0 < alpha & abs(log2FC_T0) >= lfc) |
    (padj_T1 < alpha & abs(log2FC_T1) >= lfc) |
    (padj_T2 < alpha & abs(log2FC_T2) >= lfc)
]
gene.list.DSC120 <- sig.genes.DSC120$gene
```



```{r}
#subset dds for heatmap
dds.DSC120 <- dds[, colData(dds)$geno == "DSC120_late"]

vsd.DSC120 <- vst(dds.DSC120, blind = T)
mat.DSC120 <- assay(vsd.DSC120)

keep.DSC120 <- intersect(gene.list.DSC120, rownames(mat.DSC120))
mat.DSC120  <- mat.DSC120[keep.DSC120, , drop=FALSE]

anno.col.DSC120 <- data.frame(
  Time      = colData(dds.DSC120)$time,
  Condition = colData(dds.DSC120)$condition
)
rownames(anno.col.DSC120) <- colnames(dds.DSC120)

sample.order.DSC120 <- order(colData(dds.DSC120)$time,
                             colData(dds.DSC120)$condition)

mat.DSC120 <- mat.DSC120[, sample.order.DSC120, drop=FALSE]
anno.col.DSC120 <- anno.col.DSC120[sample.order.DSC120, , drop=FALSE]

pheatmap(mat.DSC120,
         scale = "row",
         annotation_col = anno.col.DSC120,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "DSC120_late: significant genes")

```

### Clustering genes based on expression pattern - Fig 1.d 

```{r}
mat.DSC120

#scale 
mat.scaled.DSC120 <- t(scale(t(mat.DSC120)))

mat.DSC120

#scale 
mat.scaled.DSC120 <- t(scale(t(mat.DSC120)))

#check if there are NAs
any(is.na(mat.scaled.DSC120))
any(is.na(mat.DSC120))

##coldata
coldata.DSC120 <- as.data.table(colData(dds.DSC120), keep.rownames = "sample")

sample.order.DSC120 <- order(coldata.DSC120$time, coldata.DSC120$condition)

coldata.DSC120 <- coldata.DSC120[sample.order.DSC120]

#k-means clustering 
set.seed(123)
km.DSC120 <- kmeans(mat.scaled.DSC120, centers = 6)

#stable clustering
centers.DSC120 <- km.DSC120$centers
#order by mean expression 
cluster.order.DSC120 <- order(rowMeans(centers.DSC120))
#map
new.cluster.labels.DSC120 <- match(km.DSC120$cluster, cluster.order.DSC120)

#cluster dataframe
cluster.df.DSC120 <- data.table(
  gene = rownames(mat.scaled.DSC120),
  cluster = factor(new.cluster.labels.DSC120, levels = 1:6)
)

#order matrix by cluster
ord.DSC120 <- order(cluster.df.DSC120$cluster)
mat.scaled.DSC120 <- mat.scaled.DSC120[ord.DSC120, ]
cluster.df.DSC120 <- cluster.df.DSC120[ord.DSC120]

#row annotata
anno.row.DSC120 <- data.frame(Cluster = cluster.df.DSC120$cluster)
rownames(anno.row.DSC120) <- cluster.df.DSC120$gene


#heatmap
pheatmap(mat.scaled.DSC120,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anno.col.DSC120,
         annotation_row =anno.row.DSC120,
         show_rownames = F,
         main = "DSC120",
         cutree_rows = 6)

```


## DSC196_non  - Light vs NoLight - timepoints (T0,T1,T2)

```{r}
res.DSC196.T0 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC196_non")))

res.DSC196.T1 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC196_non",
                           "conditionLight.timeT1",
                           "conditionLight.genoDSC196_non.timeT1")))

res.DSC196.T2 <- results(dds,
                         contrast = list(c(
                           "condition_Light_vs_NoLight",
                           "conditionLight.genoDSC196_non",
                           "conditionLight.timeT2",
                           "conditionLight.genoDSC196_non.timeT2")))
```


```{r}
#combine
res.DSC196.all <- data.table(
  gene = rownames(res.DSC196.T0),
  
  log2FC_T0   = res.DSC196.T0$log2FoldChange,
  padj_T0     = res.DSC196.T0$padj,
  
  log2FC_T1   = res.DSC196.T1$log2FoldChange,
  padj_T1     = res.DSC196.T1$padj,
  
  log2FC_T2   = res.DSC196.T2$log2FoldChange,
  padj_T2     = res.DSC196.T2$padj
)

```


```{r}
count.T0 <- sum(res.DSC196.all$padj_T0 < alpha & abs(res.DSC196.all$log2FC_T0) >= lfc, na.rm = TRUE)
count.T1 <- sum(res.DSC196.all$padj_T1 < alpha & abs(res.DSC196.all$log2FC_T1) >= lfc, na.rm = TRUE)
count.T2 <- sum(res.DSC196.all$padj_T2 < alpha & abs(res.DSC196.all$log2FC_T2) >= lfc, na.rm = TRUE)

# 
deg.counts.DSC196 <- data.table(
  Time = c("T0", "T1", "T2"),
  N_DEGs = c(count.T0, count.T1, count.T2)
)

deg.counts.DSC196
```

```{r}
#significant DEGS
sig.genes.DSC196 <- res.DSC196.all[
  (padj_T0 < alpha & abs(log2FC_T0) >= lfc) |
    (padj_T1 < alpha & abs(log2FC_T1) >= lfc) |
    (padj_T2 < alpha & abs(log2FC_T2) >= lfc)
]
gene.list.DSC196 <- sig.genes.DSC196$gene
```



```{r}

#subset for heatmap
dds.DSC196 <- dds[, colData(dds)$geno == "DSC196_non"]

vsd.DSC196 <- vst(dds.DSC196, blind = T)
mat.DSC196 <- assay(vsd.DSC196)

keep.DSC196 <- intersect(gene.list.DSC196, rownames(mat.DSC196))
mat.DSC196  <- mat.DSC196[keep.DSC196, , drop=FALSE]

anno.col.DSC196 <- data.frame(
  Time      = colData(dds.DSC196)$time,
  Condition = colData(dds.DSC196)$condition
)
rownames(anno.col.DSC196) <- colnames(dds.DSC196)

sample.order.DSC196 <- order(colData(dds.DSC196)$time,
                             colData(dds.DSC196)$condition)

mat.DSC196 <- mat.DSC196[, sample.order.DSC196, drop=FALSE]
anno.col.DSC196 <- anno.col.DSC196[sample.order.DSC196, , drop=FALSE]

pheatmap(mat.DSC196,
         scale = "row",
         annotation_col = anno.col.DSC196,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "DSC196_non: significant genes")

```

### Clustering genes based on expression pattern - Fig 1.f
```{r}
#mat.DSC196

#scale 
mat.scaled.DSC196 <- t(scale(t(mat.DSC196)))

#scale 
mat.scaled.DSC196 <- t(scale(t(mat.DSC196)))

#check if there are NAs
any(is.na(mat.scaled.DSC196))
any(is.na(mat.DSC196))

##coldata
coldata.DSC196 <- as.data.table(colData(dds.DSC196), keep.rownames = "sample")

sample.order.DSC196 <- order(coldata.DSC196$time, coldata.DSC196$condition)

coldata.DSC196 <- coldata.DSC196[sample.order.DSC196]

#k-means clustering 
set.seed(123)
km.DSC196 <- kmeans(mat.scaled.DSC196, centers = 6)

#stable clustering
centers.DSC196 <- km.DSC196$centers
#order by mean expression 
cluster.order.DSC196 <- order(rowMeans(centers.DSC196))
#map
new.cluster.labels.DSC196 <- match(km.DSC196$cluster, cluster.order.DSC196)

#cluster dataframe
cluster.df.DSC196 <- data.table(
  gene = rownames(mat.scaled.DSC196),
  cluster = factor(new.cluster.labels.DSC196, levels = 1:6)
)

#order matrix by cluster
ord.DSC196 <- order(cluster.df.DSC196$cluster)
mat.scaled.DSC196 <- mat.scaled.DSC196[ord.DSC196, ]
cluster.df.DSC196 <- cluster.df.DSC196[ord.DSC196]

#row annotata
anno.row.DSC196 <- data.frame(Cluster = cluster.df.DSC196$cluster)
rownames(anno.row.DSC196) <- cluster.df.DSC196$gene


#heatmap
pheatmap(mat.scaled.DSC196,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anno.col.DSC196,
         annotation_row =anno.row.DSC196,
         show_rownames = F,
         main = "DSC196",
         cutree_rows = 6)

```



## Plotting the counts of signinficant DEGs

```{r}
#significan DEGS summary

degs <- rbindlist(list(
  "DSC272 (Early)" = deg.counts.DSC272,
  "DSC120 (Late)" = deg.counts.DSC120,
  "DSC196 (Non)" = deg.counts.DSC196
), idcol = "Genotype")

#long tabke
degs <- dcast(degs, Genotype ~ Time, value.var = "N_DEGs")

degs.long <- melt(degs, 
                 id.vars = "Genotype", 
                 variable.name = "Timepoint", 
                 value.name = "DEGs")

#fact
degs.long[, Genotype := factor(Genotype, levels = c("DSC272 (Early)","DSC120 (Late)","DSC196 (Non)"))]

#plot
ggplot(degs.long, aes(x = Genotype, y = DEGs, fill = Timepoint)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label =DEGs),
            position = position_dodge(width = 0.9),
            vjust = -0.1,size = 3)+
  theme_classic()+scale_fill_brewer(palette = 1) 

```



```{r}

#now maybe just to see the total counts of up and down regulated significant genes 
count296 <- res.DSC272.all[, .(
  Genotype = "DSC272 (Early)",
  Up = c(sum(padj_T0 < alpha & log2FC_T0 >= lfc, na.rm=T), 
         sum(padj_T1 < alpha & log2FC_T1 >= lfc, na.rm=T), 
         sum(padj_T2 < alpha & log2FC_T2 >= lfc, na.rm=T)),
  Down = c(sum(padj_T0 < alpha & log2FC_T0 <= -lfc, na.rm=T), 
           sum(padj_T1 < alpha & log2FC_T1 <= -lfc, na.rm=T), 
           sum(padj_T2 < alpha & log2FC_T2 <= -lfc, na.rm=T)),
  Time = c("T0", "T1", "T2")
)]

count120 <- res.DSC120.all[, .(
  Genotype = "DSC120 (Late)",
  Up = c(sum(padj_T0 < alpha & log2FC_T0 >= lfc, na.rm=T), 
         sum(padj_T1 < alpha & log2FC_T1 >= lfc, na.rm=T), 
         sum(padj_T2 < alpha & log2FC_T2 >= lfc, na.rm=T)),
  Down = c(sum(padj_T0 < alpha & log2FC_T0 <= -lfc, na.rm=T), 
           sum(padj_T1 < alpha & log2FC_T1 <= -lfc, na.rm=T), 
           sum(padj_T2 < alpha & log2FC_T2 <= -lfc, na.rm=T)),
  Time = c("T0", "T1", "T2")
)]

count196 <- res.DSC196.all[, .(
  Genotype = "DSC196 (Non)",
  Up = c(sum(padj_T0 < alpha & log2FC_T0 >= lfc, na.rm=T), 
         sum(padj_T1 < alpha & log2FC_T1 >= lfc, na.rm=T), 
         sum(padj_T2 < alpha & log2FC_T2 >= lfc, na.rm=T)),
  Down = c(sum(padj_T0 < alpha & log2FC_T0 <= -lfc, na.rm=T), 
           sum(padj_T1 < alpha & log2FC_T1 <= -lfc, na.rm=T), 
           sum(padj_T2 < alpha & log2FC_T2 <= -lfc, na.rm=T)),
  Time = c("T0", "T1", "T2")
)]

counts.stats <- rbind(count296, count120, count196)

counts.stats[, Total := Up + Down]
counts.stats <- dcast(counts.stats, Genotype ~ Time, value.var = c("Up", "Down", "Total"))
#colmun order
setcolorder(counts.stats, c("Genotype", 
                           "Up_T0", "Down_T0", "Total_T0", 
                           "Up_T1", "Down_T1", "Total_T1", 
                           "Up_T2", "Down_T2", "Total_T2"))

counts.stats
```

## Fig 1.a

```{r}
#plot counts 
#reshape
plot.dat <- melt(counts.stats, 
                  id.vars = "Genotype", 
                  measure.vars = patterns("^Up_", "^Down_"),
                  variable.name = "Timepoint",
                  value.name = c("Up", "Down"))

plot.dat[, Timepoint := factor(Timepoint, labels = c("T0", "T1", "T2"))]
#melt
plot.dat.long <- melt(plot.dat, 
                       id.vars = c("Genotype", "Timepoint"), 
                       measure.vars = c("Up", "Down"),
                       variable.name = "Exp", 
                       value.name = "Count")

#re level - i like to see early, late .. non 
plot.dat.long[, Genotype := factor(Genotype, levels = c("DSC272 (Early)","DSC120 (Late)","DSC196 (Non)"))]
#plot
ggplot(plot.dat.long, aes(x = Timepoint, y = Count, fill = Exp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label =Count),
            position = position_dodge(width = 0.9),
            vjust = -0.1,size = 2.5) +
  facet_wrap(~Genotype) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(title = "DEGs",
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

```





