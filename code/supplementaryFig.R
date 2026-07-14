library(ggplot2)
load("data/supplementary.RData")
#----SupplementaryFig2--------
FDR <- read.csv("data/DM-refference-FDR.csv", row.names = 1) %>% as.matrix()
Power <- read.csv("data/DM-refference-Power.csv", row.names = 1) %>% as.matrix()
proportion <- c("10%","20%","30%","40%","50%")
colors <- c("#7b95c6", "#49c2d9", "#67a583","#f59c7c","#c85e62")
data1 <- t(FDR)
data3 <- t(Power)
n <- ncol(data1)
m <- nrow(data1)
data2 <- data.frame(
  x = rep(c(1:4), n),
  y = as.vector(data1),
  group = rep(proportion, each=m)
)
data4 <- data.frame(
  x = rep(c(1:4), n),
  y = as.vector(data3),
  group = rep(proportion, each=m)
)

#----FDR--------
ggplot(data2, aes(x=x, y=y, color = factor(group), group=group)) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = factor(group)), size = 2, stroke = 1) +
  labs(x = "Beta", y = "FDR",title="Dirichlet-Multinomial") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 8, 2, 0)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 1),
    strip.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )


#----POWER---------
ggplot(data4, aes(x=x, y=y, color = factor(group), group=group)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = factor(group)), size = 2, stroke = 1) +
  labs(x = "Beta", y = "Reference Taxa Proportion",title="Dirichlet-Multinomial") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 8, 2, 0)) + 
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 1),
    strip.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )




#----SupplementaryFig3--------
AUPR1 <- read.csv("data/Network estimate methods comparison-AUPR.csv",row.names = 1)%>%as.matrix()
cfdr1 <- read.csv("data/Network estimate methods comparison-controlFDR.csv",row.names = 1)%>%as.matrix()
colors <- c("#c85e62","#f59c7c","#67a583","#49c2d9","#7b95c6")

data1 <- AUPR1
data2 <- data.frame(
  x = rep(c(10,20,30,40), each = 5),
  y = as.vector(data1),
  group = rep(c("BONOBO-Precision","Pearson","Spiec-Easi","SparCC","BEEM-Static"), 4)
)
target_order <- c("BONOBO-Precision","Pearson","Spiec-Easi","SparCC","BEEM-Static")
data2$group <- factor(data2$group, levels = target_order)
p1 <- ggplot(data2, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x ="Number of Taxa", y = "AUPR", title="gLV simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,17,15,18,8)) + 
  scale_y_continuous(limits = c(0.3, 0.9), breaks = seq(0.3, 0.9, 0.1)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )
p1

data3 <- cfdr1
data4 <- data.frame(
  x = rep(c(0.2,0.4,0.6,0.8), each = 5),
  y = as.vector(data3),
  group = rep(c("BONOBO-Precision","Pearson","Spiec-Easi","SparCC","BEEM-Static"), 4)
)
target_order <- c("BONOBO-Precision","Pearson","Spiec-Easi","SparCC","BEEM-Static")
data4$group <- factor(data4$group, levels = target_order)
p2 <- ggplot(data4, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x ="Recall", y = "FDR", title="gLV simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,17,15,18,8)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )
p2


methods <- c("BONOBO-Precision","SparCC","Pearson","BEEM-Static","Spiec-Easi")
data1 <- data.frame(
  Group = rep(c("BEEM-Static", "Spiec-Easi","SparCC","Pearson","BONOBO-Precision"), each = 100),
  Value = as.vector(AUC)
)
data1$Group <- factor(data1$Group, levels = methods)

ggplot(data1, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() + 
  labs( x = "Method", y = "AUC", title = "gLV") +
  theme_bw() +
  scale_fill_manual(values = c("#c85e62","#49c2d9","#f59c7c", "#7b95c6","#67a583")) +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )

methods <- c("BONOBO","SparCC","Pearson","BEEM-Static","Spiec-Easi")
data2 <- data.frame(
  Group = rep(c("BEEM-Static", "Spiec-Easi","SparCC","Pearson","BONOBO"), each = 100),
  Value = as.vector(DIFF)
)
data2$Group <- factor(data2$Group, levels = methods)

# 2. 绘图
ggplot(data2, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() + 
  labs( x = "Method", y = "RMSE", title = "gLV") +
  theme_bw() +
  scale_fill_manual(values =c("#c85e62","#49c2d9","#f59c7c", "#7b95c6","#67a583")) +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )



#----SupplementaryFig4--------
library(openxlsx)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gghalves) 
########################
# V6CRC  Zeller G(2014)
########################
df <- df1
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3
row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.V6CRC1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk1
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V6CRC2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6), 
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
    axis.title = element_text(color = "black")
  )

p.V6CRC3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.1), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Zeller G(2014) - CRC",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V6CRC3 | p.V6CRC1 | p.V6CRC2)+plot_layout(widths = c(1.5, 6, 3.5))


########################
# V19CRC  Feng Q(2015)
########################
df <- df2
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.V19CRC1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk2
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V19CRC2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6),  
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
    axis.text.y = element_blank()
  )

p.V19CRC3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.2), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Feng Q(2015) - CRC",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V19CRC3 | p.V19CRC1 | p.V19CRC2)+plot_layout(widths = c(1.5, 6, 3.5))



########################
# V43CRC  Vogtmann E(2016)
########################
df <- df3
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.V43CRC1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk3
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V43CRC2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) + 
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
    axis.title = element_text(color = "black")
  )

p.V43CRC3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.1), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Vogtmann E(2016) - CRC",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V43CRC3 | p.V43CRC1 | p.V43CRC2)+plot_layout(widths = c(1.5, 6, 3.5))



########################
# V48CD  Schirmer M(2018)
########################
df <- df4
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.V48CD1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk4
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V48CD2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)), 
    breaks = c(0,0.2,0.4,0.6),  
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
    axis.title = element_text(color = "black")
  )

p.V48CD3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.1), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Schirmer M(2018) - Crohn's Disease",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V48CD3 | p.V48CD1 | p.V48CD2)+plot_layout(widths = c(1.5, 6, 3.5))



########################
# V16CD  Nielsen HB(2014)
########################
df <- df5
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.V16CD1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk5
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V16CD2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6), 
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
    axis.title = element_text(color = "black")
  )

p.V16CD3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.2), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Nielsen HB(2014) - Crohn's Disease",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V16CD3 | p.V16CD1 | p.V16CD2)+plot_layout(widths = c(1.5, 6, 3.5))


########################
# S7CD  He Q(2017)
########################
df <- df6
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(4.5, 6.5)

p.S7CD1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )

nk <- nk6
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.S7CD2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6),  
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
    axis.title = element_text(color = "black")
  )

p.S7CD3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.3), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "He Q(2017) - Crohn's Disease",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.S7CD3 | p.S7CD1 | p.S7CD2)+plot_layout(widths = c(1.5, 6, 3.5))





########################
# V12Obesity  Liu RX(2017)
########################
df <- df7
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted2 <- melt(df, id.vars = "Method",measure.vars = "GMrepo.Hit")
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r3 <- rank(df$GMrepo.Hit)
r <- r1+r2+r3

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(3.5, 5.5)

p.V12Obesity1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  new_scale_fill() +
  
  # Hit time
  geom_tile(data = melted2,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#edb8b0","#e69191","#c25759","#C83C32"),
                       name = "Hit times") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )


nk <- nk7
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V12Obesity2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) +  
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6),  
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
    axis.title = element_text(color = "black")
  )

p.V12Obesity3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.2), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Liu RX(2017) - Obesity",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V12Obesity3 | p.V12Obesity1 | p.V12Obesity2)+plot_layout(widths = c(1.5, 6, 3.5))



########################
# V2ACVD  Jie ZY(2017)
########################
df <- df8
labelname <- colnames(df[,-1])
melted1 <- melt(df, id.vars = "Method",measure.vars = setdiff(names(df), c("Method", "GMrepo.Hit", "DA.Taxa.Prop","origin")))
melted3 <- melt(df, id.vars = "Method",measure.vars = "DA.Taxa.Prop")

r1 <- rank(df$meanTPR)
r2 <- rank(df$unionTPR)
r <- r1+r2

row_order <- df$Method[order(r)]
melted1$Method <- factor(melted1$Method, levels = row_order)
melted3$Method <- factor(melted3$Method, levels = row_order)

group_breaks <- c(3.5)

p.V2ACVD1 <- ggplot() +
  # TPR
  geom_tile(data = melted1,
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1,
            aes(x = variable, y = Method, label = round(value, 2)),
            color = "black") +
  scale_fill_gradientn(colors = c("#cce4ef","#aecfd4","#92b5ca","#599cb4","#036197"),
                       name = "TPR") +
  
  labs(x = NULL,
       y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(color = "black")
  )

nk <- nk8
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
hitdist.norm$Method <- factor(methodsname, levels = row_order)
melted <- melt(hitdist.norm, id.vars = "Method")

p.V2ACVD2 <- ggplot(melted, aes(x = variable, y = Method)) +
  geom_point(aes(size = value, color = value), alpha = 0.7) +
  scale_size(range = c(0, 6)) + 
  scale_colour_gradientn(
    colors = c("#cbead1","#a8e0b7","#73c088","#397d54","#235d3a"), 
    values = scales::rescale(c(0,0.2,0.4,0.6,1)),  
    breaks = c(0,0.2,0.4,0.6),  
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
    axis.title = element_text(color = "black")
  )

p.V2ACVD3 <- ggplot(melted3, aes(x = Method, y = -value, fill = value)) +
  geom_col() +
  geom_text(
    aes(y = -value, label = round(value, 2)),
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
    limits = c(-(max(melted3$value)+0.25), 0),
    labels = function(x) abs(x)
  )+
  theme_minimal() +
  labs(title = "Jie ZY(2017) - ACVD",
       x = NULL,
       y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  expand_limits(y = -max(melted3$value) * 1.1)

(p.V2ACVD3 | p.V2ACVD1 | p.V2ACVD2)+plot_layout(widths = c(1.5, 6, 3.5))





