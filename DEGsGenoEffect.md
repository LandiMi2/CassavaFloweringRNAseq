
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
## At time T0 

```{r}
# late vs DSC272 at T0 (NoLight)
T0.120vs272 <- results(dds, name="geno_DSC120_late_vs_DSC272_early")

# non vs early at T0 (NoLight)
T0.196vs272 <- results(dds, name="geno_DSC196_non_vs_DSC272_early")

# late vs non at T0 (NoLight)
T0.120vs196 <- results(dds, contrast=list("geno_DSC120_late_vs_DSC272_early", 
                                               "geno_DSC196_non_vs_DSC272_early"))



```



```{r}
#combine time points results 

T0.geno <- data.table(
  gene = rownames(T0.120vs272),

  log2FC_120vs272 = T0.120vs272$log2FoldChange,
  padj_120vs272   = T0.120vs272$padj,

  log2FC_196vs272 = T0.196vs272$log2FoldChange,
  padj_196vs272   = T0.196vs272$padj,

  log2FC_120vs196 = T0.120vs196$log2FoldChange,
  padj_120vs196   = T0.120vs196$padj
)

head(T0.geno)
```


```{r}

count.120vs272 <- sum(T0.geno$padj_120vs272 < alpha & abs(T0.geno$log2FC_120vs272) >= lfc, na.rm = TRUE)
count.196vs272 <- sum(T0.geno$padj_196vs272 < alpha & abs(T0.geno$log2FC_196vs272) >= lfc, na.rm = TRUE)
count.120vs196 <- sum(T0.geno$padj_120vs196 < alpha & abs(T0.geno$log2FC_120vs196) >= lfc, na.rm = TRUE)

# 
deg.counts.T0 <- data.table(
  Genos = c("120vs272", "196vs272", "120vs196"),
  N_DEGs = c(count.120vs272, count.196vs272, count.120vs196)
)

deg.counts.T0

```


```{r}
sig.genes.T0 <- T0.geno[(padj_120vs272  < alpha & abs(log2FC_120vs272 ) >= lfc) |
                                     (padj_196vs272 < alpha & abs(log2FC_196vs272) >= lfc) |
                                     (padj_120vs196 < alpha & abs(log2FC_120vs196) >= lfc)]

gene.list.T0.geno <- sig.genes.T0$gene


```


### clustering heatmap - Fig. 2b
```{r}
#any(is.na(mat.scaled.T0))
mat.T0 <- as.matrix(sig.genes.T0[, .(log2FC_120vs272, log2FC_196vs272, log2FC_120vs196)])

rownames(mat.T0) <- sig.genes.T0$gene

colnames(mat.T0) <- c("120vs272", "196vs272", "120vs196")

#scale values
mat.scaled.T0 <- t(scale(t(mat.T0)))

set.seed(123)
km.T0 <- kmeans(mat.scaled.T0, centers = 6, nstart = 25)  

#stable clustering
centers.T0 <- km.T0$centers
#order by mean expression 
cluster.order.T0 <- order(rowMeans(centers.T0))
#map
new.cluster.labels.T0 <- match(km.T0$cluster, cluster.order.T0)


#cluster dataframe
cluster.df.T0 <- data.table(
  gene = rownames(mat.scaled.T0),
  cluster = factor(new.cluster.labels.T0, levels = 1:6)
)

#order matrix by cluster
ord.T0 <- order(cluster.df.T0$cluster)
mat.scaled.T0 <- mat.scaled.T0[ord.T0, ]
cluster.df.T0 <- cluster.df.T0[ord.T0]

#row annotata
anno.row.T0 <- data.frame(Cluster = cluster.df.T0$cluster)
rownames(anno.row.T0) <- cluster.df.T0$gene

#heatmap
pheatmap(mat.scaled.T0,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = anno.row.T0,
         show_rownames = FALSE,
         main = "T0",
         angle_col = 0)


```



## At T1

```{r}
# late vs DSC272 at T1 (NoLight)
T1.120vs272 <- results(dds, contrast=list(c("geno_DSC120_late_vs_DSC272_early", "genoDSC120_late.timeT1")))

# non vs early at T1 (NoLight)
T1.196vs272 <- results(dds, contrast=list(c("geno_DSC196_non_vs_DSC272_early", "genoDSC196_non.timeT1")))

# late vs non at T1 (NoLight)
T1.120vs196 <- results(dds, contrast=list(c("geno_DSC120_late_vs_DSC272_early", "genoDSC120_late.timeT1"),
                                               c("geno_DSC196_non_vs_DSC272_early", "genoDSC196_non.timeT1")))



```



```{r}
T1.geno <- data.table(
  gene = rownames(T1.120vs272),

  log2FC_120vs272 = T1.120vs272$log2FoldChange,
  padj_120vs272   = T1.120vs272$padj,

  log2FC_196vs272 = T1.196vs272$log2FoldChange,
  padj_196vs272   = T1.196vs272$padj,

  log2FC_120vs196 = T1.120vs196$log2FoldChange,
  padj_120vs196   = T1.120vs196$padj
)

head(T1.geno)
```



```{r}
count.120vs272 <- sum(T1.geno$padj_120vs272 < alpha & abs(T1.geno$log2FC_120vs272) >= lfc, na.rm = TRUE)
count.196vs272 <- sum(T1.geno$padj_196vs272 < alpha & abs(T1.geno$log2FC_196vs272) >= lfc, na.rm = TRUE)
count.120vs196 <- sum(T1.geno$padj_120vs196 < alpha & abs(T1.geno$log2FC_120vs196) >= lfc, na.rm = TRUE)

# 
deg.counts.T1 <- data.table(
  Genos = c("120vs272", "196vs272", "120vs196"),
  N_DEGs = c(count.120vs272, count.196vs272, count.120vs196)
)

deg.counts.T1

```



```{r}
sig.genes.T1 <- T1.geno[(padj_120vs272  < alpha & abs(log2FC_120vs272 ) >= lfc) |
                                     (padj_196vs272 < alpha & abs(log2FC_196vs272) >= lfc) |
                                     (padj_120vs196 < alpha & abs(log2FC_120vs196) >= lfc)]

gene.list.T1.geno <- sig.genes.T1$gene


```

### clustering heatmap - Fig. 2c


```{r}

mat.T1 <- as.matrix(sig.genes.T1[, .(log2FC_120vs272, log2FC_196vs272, log2FC_120vs196)])

rownames(mat.T1) <- sig.genes.T1$gene

colnames(mat.T1) <- c("120vs272", "196vs272", "120vs196")

#scale values
mat.scaled.T1 <- t(scale(t(mat.T1)))

set.seed(123)
km.T1 <- kmeans(mat.scaled.T1, centers = 6, nstart = 25)  

#stable clustering
centers.T1 <- km.T1$centers
#order by mean expression 
cluster.order.T1 <- order(rowMeans(centers.T1))
#map
new.cluster.labels.T1 <- match(km.T1$cluster, cluster.order.T1)

#cluster dataframe
cluster.df.T1 <- data.table(
  gene = rownames(mat.scaled.T1),
  cluster = factor(new.cluster.labels.T1, levels = 1:6)
)

#order matrix by cluster
ord.T1 <- order(cluster.df.T1$cluster)
mat.scaled.T1 <- mat.scaled.T1[ord.T1, ]
cluster.df.T1 <- cluster.df.T1[ord.T1]

#row annotata
anno.row.T1 <- data.frame(Cluster = cluster.df.T1$cluster)
rownames(anno.row.T1) <- cluster.df.T1$gene

#heatmap
pheatmap(mat.scaled.T1,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = anno.row.T1,
         show_rownames = FALSE,
         main = "T1",
         angle_col = 0)

```



## At T2

```{r}
# late vs DSC272 at T2 (NoLight)
T2.120vs272 <- results(dds, contrast=list(c("geno_DSC120_late_vs_DSC272_early", "genoDSC120_late.timeT2")))

# non vs early at T2 (NoLight)
T2.196vs272 <- results(dds, contrast=list(c("geno_DSC196_non_vs_DSC272_early", "genoDSC196_non.timeT2")))

# late vs non at T2 (NoLight)
T2.120vs196 <- results(dds, contrast=list(c("geno_DSC120_late_vs_DSC272_early", "genoDSC120_late.timeT2"),
                                               c("geno_DSC196_non_vs_DSC272_early", "genoDSC196_non.timeT2")))


```





```{r}
#combine time points results 

T2.geno <- data.table(
  gene = rownames(T2.120vs272),

  log2FC_120vs272 = T2.120vs272$log2FoldChange,
  padj_120vs272   = T2.120vs272$padj,

  log2FC_196vs272 = T2.196vs272$log2FoldChange,
  padj_196vs272   = T2.196vs272$padj,

  log2FC_120vs196 = T2.120vs196$log2FoldChange,
  padj_120vs196   = T2.120vs196$padj
)

head(T2.geno)


```


```{r}

count.120vs272 <- sum(T2.geno$padj_120vs272 < alpha & abs(T2.geno$log2FC_120vs272) >= lfc, na.rm = TRUE)
count.196vs272 <- sum(T2.geno$padj_196vs272 < alpha & abs(T2.geno$log2FC_196vs272) >= lfc, na.rm = TRUE)
count.120vs196 <- sum(T2.geno$padj_120vs196 < alpha & abs(T2.geno$log2FC_120vs196) >= lfc, na.rm = TRUE)

# 
deg.counts.T2 <- data.table(
  Genos = c("120vs272", "196vs272", "120vs196"),
  N_DEGs = c(count.120vs272, count.196vs272, count.120vs196)
)

deg.counts.T2
```

```{r}
sig.genes.T2 <- T2.geno[(padj_120vs272  < alpha & abs(log2FC_120vs272 ) >= lfc) |
                                     (padj_196vs272 < alpha & abs(log2FC_196vs272) >= lfc) |
                                     (padj_120vs196 < alpha & abs(log2FC_120vs196) >= lfc)]

gene.list.T2.geno <- sig.genes.T2$gene


```



### clustering heatmap - Fig. 2d 
```{r}

mat.T2 <- as.matrix(sig.genes.T2[, .(log2FC_120vs272, log2FC_196vs272, log2FC_120vs196)])

rownames(mat.T2) <- sig.genes.T2$gene

colnames(mat.T2) <- c("120vs272", "196vs272", "120vs196")

#scale values
mat.scaled.T2 <- t(scale(t(mat.T2)))

set.seed(123)
km.T2 <- kmeans(mat.scaled.T2, centers = 6, nstart = 25)  

#stable clustering
centers.T2 <- km.T2$centers
#order by mean expression 
cluster.order.T2 <- order(rowMeans(centers.T2))
#map
new.cluster.labels.T2 <- match(km.T2$cluster, cluster.order.T2)

#cluster dataframe
cluster.df.T2 <- data.table(
  gene = rownames(mat.scaled.T2),
  cluster = factor(new.cluster.labels.T2, levels = 1:6)
)

#order matrix by cluster
ord.T2 <- order(cluster.df.T2$cluster)
mat.scaled.T2 <- mat.scaled.T2[ord.T2, ]
cluster.df.T2 <- cluster.df.T2[ord.T2]

#row annotata
anno.row.T2 <- data.frame(Cluster = cluster.df.T2$cluster)
rownames(anno.row.T2) <- cluster.df.T2$gene

#heatmap
pheatmap(mat.scaled.T2,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = anno.row.T2,
         show_rownames = FALSE,
         main = "T2",
         angle_col = 0)


```

A bar plot of DEGs counts

```{r}

deg.counts.T0



#significan DEGS summary
degs <- rbindlist(list(
  "T0" = deg.counts.T0,
  "T1" = deg.counts.T1,
  "T2" = deg.counts.T2
), idcol = "TimePoint")

#long tabke
degs <- dcast(degs, TimePoint ~ Genos, value.var = "N_DEGs")

degs.long <- melt(degs, 
                 id.vars = "TimePoint", 
                 variable.name = "Comparions", 
                 value.name = "DEGs")

#fact
degs.long[, TimePoint := factor(TimePoint, levels = c("T0","T1","T2"))]
degs.long[, Comparions := factor(Comparions, levels = c("120vs272","196vs272","120vs196"))]

#plot
ggplot(degs.long, aes(x = TimePoint, y = DEGs, fill = Comparions)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label =DEGs),
            position = position_dodge(width = 0.9),
            vjust = -0.1,size = 3)+
  theme_classic()+scale_fill_brewer(palette = 1)

```


Up and down regulated gene

```{r}
countT0 <- T0.geno[, .(
  TimePoint = "T0",
  Up = c(sum(padj_120vs272 < alpha & log2FC_120vs272 >= lfc, na.rm=T), 
         sum(padj_196vs272 < alpha & log2FC_196vs272 >= lfc, na.rm=T), 
         sum(padj_120vs196 < alpha & log2FC_120vs196 >= lfc, na.rm=T)),
  Down = c(sum(padj_120vs272 < alpha & log2FC_120vs272 <= -lfc, na.rm=T), 
           sum(padj_196vs272 < alpha & log2FC_196vs272 <= -lfc, na.rm=T), 
           sum(padj_120vs196 < alpha & log2FC_120vs196 <= -lfc, na.rm=T)),
  Comparison = c("120vs272", "196vs272", "120vs196")
)]

countT1 <- T1.geno[, .(
  TimePoint = "T1",
  Up = c(sum(padj_120vs272 < alpha & log2FC_120vs272 >= lfc, na.rm=T), 
         sum(padj_196vs272 < alpha & log2FC_196vs272 >= lfc, na.rm=T), 
         sum(padj_120vs196 < alpha & log2FC_120vs196 >= lfc, na.rm=T)),
  Down = c(sum(padj_120vs272 < alpha & log2FC_120vs272 <= -lfc, na.rm=T), 
           sum(padj_196vs272 < alpha & log2FC_196vs272 <= -lfc, na.rm=T), 
           sum(padj_120vs196 < alpha & log2FC_120vs196 <= -lfc, na.rm=T)),
  Comparison = c("120vs272", "196vs272", "120vs196")
)]

countT2 <- T2.geno[, .(
  TimePoint = "T2",
  Up = c(sum(padj_120vs272 < alpha & log2FC_120vs272 >= lfc, na.rm=T), 
         sum(padj_196vs272 < alpha & log2FC_196vs272 >= lfc, na.rm=T), 
         sum(padj_120vs196 < alpha & log2FC_120vs196 >= lfc, na.rm=T)),
  Down = c(sum(padj_120vs272 < alpha & log2FC_120vs272 <= -lfc, na.rm=T), 
           sum(padj_196vs272 < alpha & log2FC_196vs272 <= -lfc, na.rm=T), 
           sum(padj_120vs196 < alpha & log2FC_120vs196 <= -lfc, na.rm=T)),
  Comparison = c("120vs272", "196vs272", "120vs196")
)]


counts.stats <- rbind(countT0, countT1, countT2)

counts.stats[, Total := Up + Down]
counts.stats <- dcast(counts.stats, TimePoint ~ Comparison, value.var = c("Up", "Down", "Total"))
#colmun order
setcolorder(counts.stats, c("TimePoint", 
                            "Up_120vs272", "Down_120vs272", "Total_120vs272",
                           "Up_196vs272", "Down_196vs272", "Total_196vs272",
                           "Up_120vs196", "Down_120vs196", "Total_120vs196"))

counts.stats
```

### Bar plot up and down regulated DEGs Fig. 2a

```{r}
#plot counts 
#reshape
plot.dat <- melt(counts.stats, 
                  id.vars = "TimePoint", 
                  measure.vars = patterns("^Up_", "^Down_"),
                  variable.name = "Comparison",
                  value.name = c("Up", "Down"))

plot.dat[, Comparison := factor(Comparison, 
                                levels = 1:3, 
                                labels = c("120vs272", "196vs272", "120vs196"))]

plot.dat[, TimePoint := factor(TimePoint, labels = c("T0", "T1", "T2"))]
#melt
plot.dat.long <- melt(plot.dat, 
                       id.vars = c("TimePoint", "Comparison"), 
                       measure.vars = c("Up", "Down"),
                       variable.name = "Exp", 
                       value.name = "Count")

#re level - i like to see early, late .. non 

ggplot(plot.dat.long, aes(x = TimePoint, y = Count, fill = Exp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label =Count),
            position = position_dodge(width = 0.9),
            vjust = -0.1,size = 2.5) +
  facet_wrap(~Comparison) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(title = "DEGs",
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


```







