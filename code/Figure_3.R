library(ggnewscale)
library(ggplot2)
library(patchwork)
#-------Figure 3A----------------------------
AIC <- read.csv("data/8datasets_AIC.csv", row.names = 1)%>%as.matrix()
RMSE <- read.csv("data/8datasets_RMSE.csv", row.names = 1)%>%as.matrix()
AIC <- round(AIC,1)
RMSE <- round(RMSE,3)
INDEX <- cbind(AIC,RMSE)
dataset <- c('Zeller G(2014)','Feng Q(2015)','Vogtmann E(2016)','Schirmer M(2018)','Nielsen HB(2014)','He Q(2017)','Liu RX(2017)','Jie ZY(2017)')
row.names(INDEX) <- dataset
row.names(AIC) <- dataset
row.names(RMSE) <- dataset
colnames(INDEX) <- c("AIC0","AIC1","RMSE0","RMSE1")

#---aic-----
colnames(AIC) <- c("None","Include")
df <- as.data.frame(AIC)
df$Dataset <- rownames(df)
melted <- melt(df, id.vars = "Dataset")

p1 <- ggplot(melted, aes(x=variable, y=Dataset, fill=value)) +
  labs( x = NULL,y = "DataSets",fill = "AIC") +
  geom_tile() +
  geom_text(aes(label=value), color="black") +
  scale_fill_gradientn(colors = c("#4575B4","#91BFDB","#E0F3F8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        axis.title = element_text(size = 14, color = "black")) 

#---RMSE-----
colnames(RMSE) <- c("None","Include")
df <- as.data.frame(RMSE)
df$Dataset <- rownames(df)
melted <- melt(df, id.vars = "Dataset")

p2 <- ggplot(melted, aes(x=variable, y=Dataset, fill=value)) +
  labs(x = NULL, y = NULL, fill="RMSE") +
  geom_tile() +
  geom_text(aes(label=value), color="black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#e69191","#C83C32")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        axis.title = element_text(size = 14, color = "black")) +
  theme(axis.text.y = element_blank()) 


(p1 | p2) +
  plot_layout(guides = "keep") 


#-------Figure 3BC----------------------------
load("data/Figure3.RData")
colors <- c("darkred","#c85e62","#f47254","#f59c7c","#ffc1a6", "#a2c986","#67a583",'#a1d8e8',"#49c2d9","#7b95c6")
data1 <- t(cfdr2)

data2 <- data.frame(
  x = rep(c(0.2,0.4,0.6,0.8), 10),
  y = as.vector(data1),
  group = rep(c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA","ANCOMBC","LOCOM","ZicoSeq",
                "ADMIC-Basic","ADMIC","ADAPT"), each=4)
)
target_order <- c("ADMIC","ADAPT","MaAsLin2","ADMIC-Basic","LinDA","ALDEx2", "Wilcoxon-CLR","ANCOMBC","ZicoSeq","LOCOM")
data2$group <- factor(data2$group, levels = target_order)
p1 <- ggplot(data2, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x = "Recall", y = "FDR", title="gLV simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,1,17,2,15,0,18,4,8,3)) + 
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
p1


data3 <- t(AUPR2)
data4 <- data.frame(
  x = rep(c(10,20,30,40), 10),
  y = as.vector(data3),
  group = rep(c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA","ANCOMBC","LOCOM","ZicoSeq",
                "ADMIC-Basic","ADMIC","ADAPT"), each=4)
)
target_order <- c("ADMIC","ADAPT","MaAsLin2","ADMIC-Basic","LinDA","ALDEx2", "Wilcoxon-CLR","ANCOMBC","ZicoSeq","LOCOM")
data4$group <- factor(data4$group, levels = target_order)
p2 <- ggplot(data4, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x = "Number of Taxa", y = "AUPR", title="gLV simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,1,17,2,15,0,18,4,8,3)) + 
  scale_y_continuous(limits = c(0.4, 0.95), breaks = seq(0.4, 0.95, 0.1)) +
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

(p1 | p2) +
  plot_layout(guides = "keep") 


#-------Figure 3D----------------------------
INDEX <- read.csv("data/phage_dataset_res.csv",row.names = 1)%>%as.matrix()
INDEX <- round(INDEX,2)
df <- as.data.frame(INDEX)
df$Method <- rownames(df)
melted1 <- melt(df[,-c(3,4)], id.vars = "Method")
melted2 <- melt(df[,-c(1,2)], id.vars = "Method")

row_order <- rev(row.names(INDEX))
melted1$Method <- factor(melted1$Method, levels = row_order)
melted1$value <- round(melted1$value,3)
melted2$Method <- factor(melted2$Method, levels = row_order)
melted2$value <- round(melted2$value,3)

group_breaks <- c(2.5)

p <- ggplot() +
  # all AUPR
  geom_tile(data = melted1[1:10,],
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1[1:10,],
            aes(x = variable, y = Method, label = round(value, 3)),
            color = "black") +
  scale_fill_gradientn(colors = c("#E0F3F8","#91BFDB","#4575B4"),
                       name = "AUPR") +
  new_scale_fill() +
  
  # all FDR
  geom_tile(data = melted1[11:20,],
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted1[11:20,],
            aes(x = variable, y = Method, label = round(value, 3)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#e69191","#C83C32"),
                       name = "FDR") +
  new_scale_fill() +
  
  # resampling AUPR
  geom_tile(data = melted2[1:10,],
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2[1:10,],
            aes(x = variable, y = Method, label = round(value, 3)),
            color = "black") +
  scale_fill_gradientn(colors = c("#E0F3F8","#91BFDB","#4575B4"),
                       name = "AUPR", guide = "none") +
  new_scale_fill() +
  
  # resampling FDR
  geom_tile(data = melted2[11:20,],
            aes(x = variable, y = Method, fill = value)) +
  geom_text(data = melted2[11:20,],
            aes(x = variable, y = Method, label = round(value, 3)),
            color = "black") +
  scale_fill_gradientn(colors = c("#f5dfdb","#e69191","#C83C32"),
                       name = "FDR", guide = "none") +
  new_scale_fill() +
  
  labs(x = NULL, y = NULL) +
  geom_vline(xintercept = group_breaks, color = "white", size = 2) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        axis.title = element_text(size = 14, color = "black")) 


p

