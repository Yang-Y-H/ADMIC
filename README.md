# ADMIC

**ADMIC** is a novel method for microbiome differential abundance analysis that explicitly incorporates microbial interactions into the modeling framework. Traditional differential analysis methods typically model microbes independently, neglecting the critical impact of these microbial interactions.
To address this limitation, we established ADMIC-Basic, a high-performance and highly flexible baseline model. Building upon this foundation, ADMIC integrates microbial interaction terms as confounders. By incorporating interaction terms, ADMIC filters out secondary differential microbes whose abundance differences are driven by microbial interactions, thereby enhancing the capability to identify primary differential microbes directly associated with the disease.


```r
# To use ADMIC, you need to load the following R packages
install.packages("stats", "glasso")
library(stats)
library(glasso)
```

## 📊 Example Data

We provide a toy **case-control dataset** for you to test the package. You can find and download the files directly from the `data/` directory:

* 📊 [`data/count_table.csv`](./data/data_count.csv) — Microbial abundance count table (rows as samples, columns as OTUs/taxa).
* 📋 [`data/metadata.csv`](./data/data_meta.csv) — Sample metadata (including case-control status and other covariates).

```r
source("ADMIC_fun.R")
source("bonobo.precision_fun.R")
source("reference.select_fun.R")
# -----------------------------------
#  Loading data
# -----------------------------------
data <- read.csv("data/data_count.csv", row.names = 1)
data <- as.matrix(data)
meta <- read.csv("data/data_meta.csv", row.names = 1)
meta$Phenotype <- factor(meta$Phenotype)
meta$Gender <- factor(meta$Gender)
meta$Population <- factor(meta$Population)

# -----------------------------------
#  Filtering taxa
# -----------------------------------
group1 <- which(meta$Phenotype == levels(meta$Phenotype)[1])
group2 <- which(meta$Phenotype == levels(meta$Phenotype)[2])
#group1
prop.presence1 <- colMeans(data[group1,]>0)
otus.keep1 <- which(prop.presence1 >=0.05)
#group2
prop.presence2 <- colMeans(data[group2,]>0)
otus.keep2 <- which(prop.presence2 >=0.05)
otus.keep <- intersect(otus.keep1, otus.keep2)
data.filter <- data[,otus.keep]

Y <- meta$Phenotype
C <- meta[, colnames(meta) != "Phenotype"]

RES <- ADMIC(data.filter,Y,C)
RES.noC <- ADMIC(data.filter,Y,C = NULL)
```
