library(openxlsx)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
#-----------Figure 7-------------------------------
data <- read.xlsx("data/Index.xlsx")
my_order <- rev(data$Method)

melted1 <- melt(data[,c(1,2)], id.vars = "Method")
melted1$Method <- factor(melted1$Method, levels = my_order)
colnames(melted1) <- c("Method","Label","TPR")

melted2 <- melt(data[,c(1,3)], id.vars = "Method")
melted2$Method <- factor(melted2$Method, levels = my_order)
colnames(melted2) <- c("Method","Label","meanTPR")

melted3 <- melt(data[,c(1,4)], id.vars = "Method")
melted3$Method <- factor(melted3$Method, levels = my_order)
colnames(melted3) <- c("Method","Label","GMrepo.Hit")

melted4 <- melt(data[,c(1,5)], id.vars = "Method")
melted4$Method <- factor(melted4$Method, levels = my_order)
colnames(melted4) <- c("Method","Label","TaxaProp")


p.heatmap1 <- ggplot() +
  # union TPR
  geom_tile(data = melted1,
            aes(x = Label, y = Method, fill = TPR)) +
  geom_text(data = melted1,
            aes(x = Label, y = Method, label = round(TPR, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#FFE5CC","#FBC993","#FCAE60","#FEA344","#FE8A16"),
                       name = "unionTPR") +
  new_scale_fill() +
  
  # mean TPR
  geom_tile(data = melted2,
            aes(x = Label, y = Method, fill = meanTPR)) +
  geom_text(data = melted2,
            aes(x = Label, y = Method, label = round(meanTPR, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),                     
                       name = "meanTPR") +
  new_scale_fill() +
  
  # GMrepo.Hit
  geom_tile(data = melted3,
            aes(x = Label, y = Method, fill = GMrepo.Hit)) +
  geom_text(data = melted3,
            aes(x = Label, y = Method, label = round(GMrepo.Hit, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),                     
                       name = "GMrepo Hit") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = 2.5, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(), 
    axis.title = element_text(color = "black"),
    legend.position = "left",     
    legend.box = "horizontal",      
    legend.box.just = "top",       
    legend.spacing.x = unit(0.25, "cm"), 
    legend.title = element_text(size = 9) 
  )

nk1 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V6CRC")
nk2 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V19CRC")
nk3 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V43CRC")
nk4 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V48CD")
nk5 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V16CD")
nk6 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "S7CD")
nk7 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V12Obesity")
nk8 <- read.xlsx("data/Figure7-plotdata.xlsx",sheet = "V2ACVD")
nk <- rbind(nk1,nk2,nk3,nk4,nk5,nk6,nk7,nk8)
nk <- as.matrix(nk)
data <- nk[,-1]
methodsname <- colnames(data)
n <- length(methodsname)
hit <- rowSums(data)
detectnum <- colSums(data)

hitdist <- matrix(0, n, n)
meanhit <- c()
for(i in 1:n){
  dataxa <- which(data[,i]==1)
  dahit <- hit[dataxa]
  meanhit[i] <- mean(dahit)
  for(j in 1:n){
    hitdist[i,j] <- length(which(dahit==j))
  }
}
hitdist.norm <- hitdist/detectnum
row.names(hitdist.norm) <- methodsname
colnames(hitdist.norm) <- c(1:n)
hitdist.norm <- as.data.frame(hitdist.norm)
hitdist.norm$Method <- factor(methodsname, levels = my_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.bubble <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
 
  scale_size(range = c(0, 8)) +  # 气泡大小范围
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    guide = "colorbar"
  ) +
  theme_minimal() +
  labs(
    x = "Hit Number",
    y = NULL,
    size = "Taxa Proportion",
    color = "Taxa Proportion"
  )+
  theme(
    axis.text.x = element_text(hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),  
    axis.title = element_text(color = "black"),
    legend.title = element_text(size = 9)
  )
p.bubble

p.bar <- ggplot(melted4, aes(x = Method, y = -TaxaProp, fill = TaxaProp)) +
  geom_col() +
  geom_text(
    aes(y = -TaxaProp, label = round(TaxaProp, 2)),
    hjust = 1.1,   
    size = 3,
    color = "black"
  ) +
  coord_flip() +
  scale_fill_gradientn(
    colors = c("#D3D3D3", "#BEBEBE", "#A9A9A9", "#696969", "#808080"),
    name = "DA Taxa Prop"
  ) +
  scale_y_continuous(
    limits = c(-0.5, 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted4$TaxaProp) * 1.1)

(p.bar | p.heatmap1 | p.bubble)+plot_layout(widths = c(3, 5, 6))


 




