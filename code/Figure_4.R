library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyr)
library(dplyr)
library(rstatix)
library(patchwork)
library(gghalves) 
load("data/Figure4.RData")
#-----------Figure 4AB-------------------------------
data <- read.csv("data/Figure4AB.csv", row.names = 1)
df1 <- data.frame(
  value = c(data$n0, data$n.eff),
  group = factor(rep(c("ADMIC-Basic","ADMIC"), each = nrow(data)))
)
df1$group <- factor(df1$group)
my_comparisons1 <- list(c("ADMIC-Basic","ADMIC")) 
pA1 <- ggplot(df1, aes(x = group, y = value, fill = group)) +
  ylim(0, max(df1$value)+2) +
  geom_boxplot(width = 0.5, alpha = 0.9,  median.colour = NA ) +
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,fatten = 1,color = "black",linewidth = 1.2) +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = "Number of differential taxa",title = "gLV simulated data") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black")) +
  scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  scale_color_manual(values = c("#F8766D","#00BFC4"))+
  stat_compare_means(
    comparisons = my_comparisons1, 
    method = "wilcox.test", 
    method.args = list(alternative = "greater"),
    label = "p.signif",          
    bracket.size = 0.3,          
    tip.length = 0.02           
  )


df2 <- data.frame(
  value = c(data$fdr0, data$fdr.adj),
  group = factor(rep(c("all detected","adjust detected"), each = nrow(data)))
)
df2$group <- factor(df2$group)
my_comparisons2 <- list(c("all detected","adjust detected")) 
pA2 <- ggplot(df2, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.9,  median.colour = NA ) + 
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,fatten = 1,color = "black",linewidth = 1.2) +
  theme_bw(base_size = 14) +
  labs(title = "gLV",x = NULL, y = "FDR") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black")) +
  scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  scale_color_manual(values = c("#F8766D","#00BFC4"))+
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1.1, by = 0.2)
  ) +
  stat_compare_means(
    comparisons = my_comparisons2, 
    method = "wilcox.test", 
    method.args = list(alternative = "less"),
    label = "p.signif",           
    bracket.size = 0.3,           
    tip.length = 0.02            
  )


#-----------Figure 4C-------------------------------
library(ggplot2)
library(dplyr)
library(rstatix)
data <- rbind(data.sim, data.real.all, data.phage)
data$ratio <- as.numeric(data$ratio)
data$ratio.log <- log(data$ratio)
DATASET <- c("gLV",'Schirmer M(2018)','Vogtmann E(2016)',"He Q(2017)","Liu RX(2017)","Zeller G(2014)",
             "Jie ZY(2017)", "Nielsen HB(2014)", 'Feng Q(2015)',"T4+F1")
n <- length(DATASET)
meanvalue <- c()
for(i in 1:n){
  meanvalue[i] <- mean(data$ratio[which(data$datasets == DATASET[i])])
}
DATASET[order(meanvalue)]
DATASET <- DATASET[order(meanvalue)]
data$datasets <- factor(data$datasets, levels =DATASET)


data1 <- data[-which(data$datasets == "gLV"),]
data1 <- data1[-which(data1$datasets == "T4+F1"),]
DATASET <- setdiff(DATASET,c("gLV","T4+F1"))
n <- length(DATASET)
maxvalue <- c()
for(i in 1:n){
  maxvalue[i] <- quantile(data1$ratio[which(data1$datasets == DATASET[i])],0.95)
}
names(maxvalue) <- DATASET
stat.test_data1 <- data1 %>%
  group_by(datasets) %>%
  filter(n() > 1 && sd(ratio, na.rm = TRUE) != 0) %>%
  t_test(ratio ~ 1, mu = 1, alternative = "greater") %>%
  add_significance("p") %>%                             
  add_xy_position(x = "Dataset", dodge = 0.8) %>%
  mutate(y.position = maxvalue+0.02, x = datasets)  

pB1 <- ggplot(data1, aes(x = datasets, y = ratio, fill = datasets)) +
  stat_summary(
    fun.data = function(x) {
      data.frame(ymin = quantile(x, probs = 0.05, na.rm = TRUE),
                 ymax = quantile(x, probs = 0.95, na.rm = TRUE))
    },
    geom = "errorbar",
    width = 0.3,
    linewidth = 1.0,
    color = "#00BFC4"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 17,        
    size = 3,
    color = "black"
  ) +
  geom_text(
    data = stat.test_data1,
    aes(x = datasets, y = y.position, label = p.signif), 
    inherit.aes = FALSE,
    size = 4, fontface = "bold",
    angle = 270,        
    hjust = 0.5,         
    vjust = 0.5 
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(
    limits = c(0.8, max(maxvalue + 0.02))
  ) +
  theme_bw(base_size = 14) +
  coord_flip() +
  labs(x = NULL, y = "Correlation coefficients Ratio") +
  theme(legend.position = "none") 


DATASET <- c("T4+F1","gLV")
data2 <- data[which(data$datasets %in% DATASET),]
data2$datasets <- factor(data2$datasets, levels =DATASET)
n <- length(DATASET)
maxvalue <- c()
for(i in 1:n){
  maxvalue[i] <- quantile(data2$ratio[which(data2$datasets == DATASET[i])],0.95)
}

stat.test_data2 <- data2 %>%
  group_by(datasets) %>%
  t_test(ratio ~ 1, mu = 1, alternative = "greater") %>%
  add_significance("p") %>%                             
  add_xy_position(x = "Dataset", dodge = 0.8) %>%
  mutate(y.position = maxvalue+0.2, x = datasets)  


pB2 <- ggplot(data2, aes(x = datasets, y = ratio, fill = datasets)) +
  stat_summary(
    aes(color = datasets),
    fun.data = function(x) {
      data.frame(ymin = quantile(x, probs = 0.05, na.rm = TRUE),
                 ymax = quantile(x, probs = 0.95, na.rm = TRUE))
    },
    geom = "errorbar",
    width = 0.3,
    linewidth = 1.0
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 17,       
    size = 3,
    color = "black"
  ) +
  geom_text(
    data = stat.test_data2,
    aes(x = datasets, y = y.position, label = p.signif), 
    inherit.aes = FALSE,
    size = 4, fontface = "bold",
    angle = 270,         
    hjust = 0.5,           
    vjust = 0.5 
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  scale_color_manual(values = c("#00BFC4","#F8766D"))+
  theme_bw(base_size = 14) +
  coord_flip() +
  labs(x =NULL, y =NULL) +
  theme(legend.position = "none") 



pB <- (pB2/pB1) +
  plot_layout(
    heights = c(2,8),
    guides = "keep"
  ) 
library(patchwork)
(pA1 | pA2 | pB) +
  plot_layout(
    widths = c(0.9,0.9,1),
    guides = "keep"
  ) 



