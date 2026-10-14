

```{r setup, include=F}
library(data.table)
library(pathview)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(org.At.tair.db)
library(magick)

#load significan DEGS
sig.272 <- fread("../sigDEGs/sig.genes.DSC272.txt")
sig.120 <- fread("../sigDEGs/sig.genes.DSC120.txt")
sig.196 <- fread("../sigDEGs/sig.genes.DSC196.txt")

#background 
universe <- fread("../DEGs/countsFilterdwithArab.txt")

```

## Run GO and KEGG pathwas on 272 genotypes (light effect) - Fig. 1c
```{r}
#run enrichGo
go.272 <- enrichGO(
  gene = sig.272$arabID,
  OrgDb = org.At.tair.db,
  universe = universe$arabID,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)
#check all terms 
go.272 <- enrichGO(
  gene = sig.272$arabID,
  OrgDb = org.At.tair.db,
  universe = universe$arabID,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr",pvalueCutoff = 1, qvalueCutoff = 1,
)

go.272.def <- enrichGO(
  gene = sig.272$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.272.dt <- as.data.table(go.272)
go.272.def.dt <- as.data.table(go.272.def)
#write.table(go.272.dt ,"../sigDEGs/GO/go.272.BG.txt",sep="\t", row.names = F, quote = F)
#write.table(go.272.def.dt ,"../sigDEGs/GO/go.272.txt",sep="\t", row.names = F, quote = F)
#write.table(go.272.dt ,"../sigDEGs/GO/go.272.ALL.txt",sep="\t", row.names = F, quote = F)

## plot only the universe one
a <- dotplot(go.272, showCategory = 20, x = "p.adjust",label_format=50,
             title = "272") 

a$data$ONTOLOGY <- factor(a$data$ONTOLOGY)

a + aes(color = ONTOLOGY)

```



## Run GO and KEGG pathwas on 120 genotypes (light effect) - Fig. 1e
```{r}
#run enrichGo
go.120 <- enrichGO(
  gene = sig.120$arabID,
  OrgDb = org.At.tair.db,
  universe = universe$arabID,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

go.120.def <- enrichGO(
  gene = sig.120$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.120.dt <- as.data.table(go.120)
go.120.def.dt <- as.data.table(go.120.def)
write.table(go.120.dt ,"../sigDEGs/GO/go.120.BG.txt",sep="\t", row.names = F, quote = F)
write.table(go.120.def.dt ,"../sigDEGs/GO/go.120.txt",sep="\t", row.names = F, quote = F)

b <- dotplot(go.120, showCategory = 20, x = "p.adjust",label_format=50,
             title = "120") 

b$data$ONTOLOGY <- factor(b$data$ONTOLOGY)

b + aes(color = ONTOLOGY)
```



## Run GO and KEGG pathwas on 196 genotypes (light effect) - Fig. 1j
```{r}
#run enrichGo
go.196 <- enrichGO(
  gene = sig.196$arabID,
  universe = universe$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

go.196.def <- enrichGO(
  gene = sig.196$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.196.dt <- as.data.table(go.196)  ## 115 Go terms 
go.196.def.dt <- as.data.table(go.196.def)
write.table(go.196.dt ,"../sigDEGs/GO/go.196.BG.txt",sep="\t", row.names = F, quote = F)
write.table(go.196.def.dt ,"../sigDEGs/GO/go.196.txt",sep="\t", row.names = F, quote = F)

c <- dotplot(go.196, showCategory = 20, x = "p.adjust",label_format=50,
             title = "196") 

c$data$ONTOLOGY <- factor(c$data$ONTOLOGY)

c + aes(color = ONTOLOGY)

```



## Now check out Genotype differences at No light

```{r}
#load 
sig.geno.T0 <- fread("../sigDEGs/geno/sig.genes.T0.exp.txt")
sig.geno.T1 <- fread("../sigDEGs/geno/sig.genes.T1.exp.txt")
sig.geno.T2 <- fread("../sigDEGs/geno/sig.genes.T2.exp.txt")


```


```{r}
#run enrichGo
go.T0 <- enrichGO(
  gene = sig.geno.T0$arabID,
  universe = universe$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

go.T0.def <- enrichGO(
  gene = sig.geno.T0$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)


#get results
go.T0.dt <- as.data.table(go.T0)
go.T0.def.dt <- as.data.table(go.T0.def)
write.table(go.T0.dt,"../sigDEGs/geno/GO/go.T0.BG.txt",sep="\t", row.names = F, quote = F)
write.table(go.T0.def.dt,"../sigDEGs/geno/GO/go.T0.txt",sep="\t", row.names = F, quote = F)

a1 <- dotplot(go.T0, showCategory = 20, x = "p.adjust",label_format=50,
             title = "T0") 

a1$data$ONTOLOGY <- factor(a1$data$ONTOLOGY)

a1 + aes(color = ONTOLOGY)


```


```{r}
go.T1 <- enrichGO(
  gene = sig.geno.T1$arabID,
  OrgDb = org.At.tair.db,
  universe = universe$arabID,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

go.T1.def <- enrichGO(
  gene = sig.geno.T1$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.T1.dt <- as.data.table(go.T1)
go.T1.def.dt <- as.data.table(go.T1.def)
write.table(go.T1.dt,"../sigDEGs/geno/GO/go.T1.BG.txt",sep="\t", row.names = F, quote = F)
write.table(go.T1.def.dt,"../sigDEGs/geno/GO/go.T1.txt",sep="\t", row.names = F, quote = F)

a2 <- dotplot(go.T1, showCategory = 20, x = "p.adjust",label_format=50,
             title = "T1") 

a2$data$ONTOLOGY <- factor(a2$data$ONTOLOGY)

a2 + aes(color = ONTOLOGY)

```


```{r}
#run enrichGo
go.T2 <- enrichGO(
  gene = sig.geno.T2$arabID,
  universe = universe$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

go.T2.def <- enrichGO(
  gene = sig.geno.T2$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.T2.dt <- as.data.table(go.T2) #38 GO terms
go.T2.def.dt <- as.data.table(go.T2.def) 
write.table(go.T2.dt,"../sigDEGs/geno/GO/go.T2.BG.txt",sep="\t", row.names = F, quote = F)
write.table(go.T2.def.dt,"../sigDEGs/geno/GO/go.T2.txt",sep="\t", row.names = F, quote = F)

a3 <- dotplot(go.T2, showCategory = 20, x = "p.adjust",label_format=50,
             title = "T2") 

a3$data$ONTOLOGY <- factor(a3$data$ONTOLOGY)

a3 + aes(color = ONTOLOGY)
```
