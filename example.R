rm(list = ls())
library(stats)
library(glasso)

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
#  Filtering otus
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




