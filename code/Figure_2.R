library(magrittr)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(reshape2)
library(patchwork)
#-------Figure 2----------------------------
#----DM simulation resluts------
AUPR <- read.csv("data/DM-Basic-AUPR.csv", row.names = 1)
FDR <- read.csv("data/DM-Basic-FDR.csv", row.names = 1)
SEN <- read.csv("data/DM-Basic-SEN.csv", row.names = 1)

colors <- c("#c85e62","#f47254","#f59c7c","#ffc1a6","#a2c986","#67a583",'#a1d8e8',"#49c2d9","#7b95c6")

#---AUPR----
data1 <- t(AUPR)
data2 <- data.frame(
  x = rep(seq(0,4,1), 9),
  y = as.vector(data1),
  group = rep(c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA","ANCOMBC","LOCOM","ZicoSeq","ADAPT",
                "ADMIC-Basic"), each=5)
)
target_order <- c("ADMIC-Basic","LOCOM","ZicoSeq","ADAPT","ALDEx2","ANCOMBC", "LinDA","MaAsLin2","Wilcoxon-CLR")
data2$group <- factor(data2$group, levels = target_order)
p1 <- ggplot(data2, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  # 4. 标签
  labs(x = "Beta", y = "AUPR", title="DM simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,1,17,2,15,0,18,4,8)) + 
  scale_y_continuous(limits = c(0.2, 0.9), breaks = seq(0.2, 0.9, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )
p1

#----FDR--------
data3 <- t(FDR)
data4 <- data.frame(
  x = rep(seq(0,4,1), 9),
  y = as.vector(data3),
  group = rep(c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA","ANCOMBC","LOCOM","ZicoSeq","ADAPT",
                "ADMIC-Basic"), each=5)
)
target_order <- c("ADMIC-Basic","LOCOM","ZicoSeq","ADAPT","ALDEx2","ANCOMBC", "LinDA","MaAsLin2","Wilcoxon-CLR")
data4$group <- factor(data4$group, levels = target_order)
p2 <- ggplot(data4, aes(x=x, y=y, color = group, group=group)) + 
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray2", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x = "Beta", y = "FDR", title="DM simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,1,17,2,15,0,18,4,8)) + 
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
(p1 | p2) +
  plot_layout(guides = "keep") 






