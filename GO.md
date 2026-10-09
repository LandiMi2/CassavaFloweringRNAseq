

```{r setup, include=F}
library(data.table)
library(pathview)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(org.At.tair.db)
library(magick)

setwd("~/Documents/PhD-Bioinformatics-EpiCass/Analysis/Flowering")
#load significan DEGS
sig.272 <- fread("sigDEGs/sig.genes.DSC272.txt")
sig.120 <- fread("sigDEGs/sig.genes.DSC120.txt")
sig.196 <- fread("sigDEGs/sig.genes.DSC196.txt")

```

## Run GO on 272 genotypes (light effect) - Fig. 1c
```{r}
#run enrichGo
go.272 <- enrichGO(
  gene = sig.272$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.272.dt <- as.data.table(go.272)
#write.table(go.272.dt ,"sigDEGs/GO/go.272.txt",sep="\t", row.names = F, quote = F)
#save 

a <- dotplot(go.272, showCategory = 20, x = "p.adjust",label_format=50,
             title = "272") 

a$data$ONTOLOGY <- factor(a$data$ONTOLOGY)

a + aes(color = ONTOLOGY)
```

## Run GO on 120 genotypes (light effect) - Fig. 1e
```{r}
#run enrichGo
go.120 <- enrichGO(
  gene = sig.120$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.120.dt <- as.data.table(go.120)
#write.table(go.120.dt ,"sigDEGs/GO/go.120.txt",sep="\t", row.names = F, quote = F)

b <- dotplot(go.120, showCategory = 20, x = "p.adjust",label_format=50,
             title = "120") 

b$data$ONTOLOGY <- factor(b$data$ONTOLOGY)

b + aes(color = ONTOLOGY)
```


## Run GO on 196 genotypes (light effect) - Fig. 1j
```{r}
#run enrichGo
go.196 <- enrichGO(
  gene = sig.196$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.196.dt <- as.data.table(go.196)  ## 115 Go terms 
#write.table(go.196.dt ,"sigDEGs/GO/go.196.txt",sep="\t", row.names = F, quote = F)

c <- dotplot(go.196, showCategory = 20, x = "p.adjust",label_format=50,
             title = "196") 

c$data$ONTOLOGY <- factor(c$data$ONTOLOGY)

c + aes(color = ONTOLOGY)

```


## Now check out Genotype differences at No light

```{r}
#load 
sig.geno.T0 <- fread("sigDEGs/geno/sig.genes.T0.exp.txt")
sig.geno.T1 <- fread("sigDEGs/geno/sig.genes.T1.exp.txt")
sig.geno.T2 <- fread("sigDEGs/geno/sig.genes.T2.exp.txt")


```

## GO and on T0 - Fig. 2e
```{r}
#run enrichGo
go.T0 <- enrichGO(
  gene = sig.geno.T0$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)


#get results
go.T0.dt <- as.data.table(go.T0)
#write.table(go.T0.dt,"sigDEGs/geno/GO/go.T0.txt",sep="\t", row.names = F, quote = F)

a1 <- dotplot(go.T0, showCategory = 20, x = "p.adjust",label_format=50,
             title = "T0") 

a1$data$ONTOLOGY <- factor(a1$data$ONTOLOGY)

a1 + aes(color = ONTOLOGY)


```


```{r}
go.T1 <- enrichGO(
  gene = sig.geno.T1$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)
# no significant enrichnemt 

```
## GO and on T2 - Fig. 2f

```{r}
#run enrichGo
go.T2 <- enrichGO(
  gene = sig.geno.T2$arabID,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "ALL",
  pAdjustMethod = "fdr"
)

#get results
go.T2.dt <- as.data.table(go.T2) #38 GO terms
#write.table(go.T2.dt,"sigDEGs/geno/GO/go.T2.txt",sep="\t", row.names = F, quote = F)

a2 <- dotplot(go.T2, showCategory = 20, x = "p.adjust",label_format=50,
             title = "T2") 

a2$data$ONTOLOGY <- factor(a2$data$ONTOLOGY)

a2 + aes(color = ONTOLOGY)
```
