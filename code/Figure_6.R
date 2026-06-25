#-----------Figure 6A-------------------------------
Netchangepor1 <- read.csv("data/Scope of Network Changes.csv",row.names = 1)%>%as.matrix()

#----Absolute difference >0.2------
df <- data.frame(
  value = c(Netchangepor1[,1], Netchangepor1[,2], Netchangepor1[,3]),
  group = factor(rep(c("SpiecEasi","SparCC","Pearson"), each = 8))
)
df$group <- factor(df$group, levels = c("SpiecEasi","SparCC","Pearson"))

p.netchange1 <- ggplot(df, aes(x = group, y = value, fill = group)) +
  ylim(-0.01, 0.7) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) + 
  theme_bw() +
  labs(x = "Methods", y = "Network Change Proportion",title = "Absolute Difference") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#599cb4","#c25759", "#FEA344")) +
  scale_color_manual(values = c("#599cb4","#c25759", "#FEA344"))
p.netchange1

#----Relative difference >0.2------
df <- data.frame(
  value = c(Netchangepor1[,4], Netchangepor1[,5], Netchangepor1[,6]),
  group = factor(rep(c("SpiecEasi","SparCC","Pearson"), each = 8))
)
df$group <- factor(df$group, levels = c("SpiecEasi","SparCC","Pearson"))

p.netchange2 <- ggplot(df, aes(x = group, y = value, fill = group)) +
  ylim(-0.01,0.7)+
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) + # 箱线图
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) + # 散点
  theme_bw() +
  labs(x = "Methods", y = "Network Change Proportion",title = "Ratio Difference")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#599cb4","#c25759", "#FEA344")) +
  scale_color_manual(values = c("#599cb4","#c25759", "#FEA344"))
p.netchange2



#-----------Figure 6B-------------------------------
load("data/Figure6.RData")

data1 <- t(AUPR2)
data2 <- data.frame(
  x = rep(seq(0,0.6,0.1), 10),
  y = as.vector(data3),
  group = rep(c("Wilcoxon-CLR","MaAsLin2", "ALDEx2","LinDA","ANCOMBC","LOCOM","ZicoSeq",
                "ADMIC-Basic","ADMIC","ADAPT"), each=7)
)
target_order <- c("ADMIC","ADAPT","MaAsLin2","ADMIC-Basic","LinDA","ALDEx2", "Wilcoxon-CLR","ANCOMBC","ZicoSeq","LOCOM")
data2$group <- factor(data2$group, levels = target_order)
p1 <- ggplot(data2, aes(x=x, y=y, color = group, group=group)) + 
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = group), size = 2, stroke = 1) +
  labs(x = "Proportion of Network Changes", y = "AUPR", title="gLV simulated data") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16,1,17,2,15,0,18,4,8,3)) + 
  scale_y_continuous(limits = c(0.5, 0.85), breaks = seq(0.5, 0.85, 0.1)) +
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

library(patchwork)
(p.netchange1 | p.netchange2 | p1) +
  plot_layout(guides = "keep") 
