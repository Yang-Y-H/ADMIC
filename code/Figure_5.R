library(tidyr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

data <- read.csv("data/Figure5B.csv",row.names = 1)
taxa.num.mat <- read.csv("data/Figure5A_taxa_num_mat.csv",row.names = 1)
data.mean <- taxa.num.mat
n <- nrow(data.mean)
for(i in 1:n){
  data.mean[i,-1] <- data[which(data$Dataset == data.mean$Dataset[i]),3:5] %>% colMeans()
}
df <- as.data.frame(data)
df$ADMIC.Basic <- as.numeric(df$ADMIC.Basic)
df$ADMIC <- as.numeric(df$ADMIC)
df$ADAPT <- as.numeric(df$ADAPT)
my_order <- c( 'Schirmer M(2018)', 'Feng Q(2015)', 
              'Zeller G(2014)','Vogtmann E(2016)',  'Jie ZY(2017)', 
              'Liu RX(2017)', 'Nielsen HB(2014)', 'He Q(2017)',
              'B40-8+VD13','T4+F1')

df$Dataset <- factor(df$Dataset, levels = my_order)

taxa.num.mat$Dataset <- factor(taxa.num.mat$Dataset, levels = my_order)
data.mean$Dataset <- factor(data.mean$Dataset, levels = my_order)

#-----------Figure 5B-------------------------------
data_long <- df %>%
  pivot_longer(
    cols = c("ADMIC.Basic", "ADMIC","ADAPT"), 
    names_to = "Metric",                   
    values_to = "Value"                    
  )
data_long$Metric <- factor(data_long$Metric , levels = c("ADAPT", "ADMIC.Basic", "ADMIC"))
data_plot <- data_long %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("ADAPT", "ADMIC.Basic", "ADMIC")
    )
  )


data_long2 <- data.mean %>%
  pivot_longer(
    cols = c("ADMIC.Basic", "ADMIC", "ADAPT"),
    names_to = "Metric",                  
    values_to = "Value"                  
  )
data_long2$Metric <- factor(data_long2$Metric , levels = c("ADAPT", "ADMIC.Basic", "ADMIC"))
data_plot2 <- data_long2 %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("ADAPT", "ADMIC.Basic", "ADMIC")
    )
  )

pairwise_res <- data_long %>%
  pairwise_wilcox_test(
    Value ~ Metric,
    paired = F,
    p.adjust.method = "none"
  ) %>%
  filter(p.adj < 0.05) 
pairwise_res <- pairwise_res %>%
  mutate(
    y.position = c(1.02, 1.05, 1.08)[row_number(pairwise_res)]
  )

p1 <- ggplot(data_plot, aes(x = Metric, y = Value)) +
  geom_boxplot(
    fill = "white",       
    colour = "black",      
    width = 0.55,
    outlier.shape = NA,
    alpha = 0,
    color = "black",
    median.colour = NA 
  ) +
  stat_summary(fun = mean,geom = "crossbar", width = 0.55,fatten = 1,color = "black",linewidth = 1) +

  geom_point(
    data = data_plot2,
    aes(x = Metric, y = Value, color = Dataset, shape = Dataset),
    size = 2,
    alpha = 1
  ) +
  scale_shape_manual(
    values = c(
      "Schirmer M(2018)" = 2,
      "Feng Q(2015)" = 3,
      "Zeller G(2014)" = 4,
      "Vogtmann E(2016)" = 5,
      "Jie ZY(2017)" = 15,
      "Liu RX(2017)" = 16,
      "Nielsen HB(2014)" = 17,
      "He Q(2017)" = 18,
      "B40-8+VD13" = 0,
      "T4+F1" = 8
    )
  ) +
  geom_line(
    data = data_plot2,
    aes(x = Metric, y = Value, group = Dataset, color = Dataset),
    linewidth = 0.6,
    alpha = 1
  )+
  stat_pvalue_manual(
    pairwise_res,
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  # y轴
  scale_y_continuous(
    limits = c(0.35, 1.08),
    breaks = seq(0.35, 1.08, by = 0.1)
  ) +
  scale_fill_manual(
    values = c(
      "ADAPT" = "#7b95c6",
      "ADMIC.Basic" = "#DB9D5F",
      "ADMIC" = "#c85e62"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "AUC",
    fill = "Metric",
    color = "Dataset"
  ) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "right"
  )
p1
  


#-----------Figure 5A-------------------------------
data_long1 <- taxa.num.mat %>%
  pivot_longer(
    cols = c("ADMIC.Basic", "ADMIC", "ADAPT"), 
    names_to = "Metric",                  
    values_to = "Value"                   
  )
data_long1$Metric <- factor(data_long1$Metric , levels = c("ADAPT", "ADMIC.Basic", "ADMIC"))

data_plot1 <- data_long1 %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("ADAPT", "ADMIC.Basic", "ADMIC")
    )
  )
pairwise_res1 <- data_long1 %>%
  pairwise_wilcox_test(
    Value ~ Metric,
    paired = T,
    p.adjust.method = "none"
  ) %>%
  filter(p.adj < 0.05) 
pairwise_res1 <- pairwise_res1%>%
  mutate(
    y.position = c(72, 76, 80)[row_number(pairwise_res1)] 
  )


p2 <- ggplot(data_plot1, aes(x = Metric, y = Value)) +
  geom_boxplot(
    fill = "white",       
    color = "black",       
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.9,
    color = "black",
    median.colour = NA 
  ) +
  stat_summary(fun = mean,geom = "crossbar", width = 0.55,fatten = 1,color = "black",linewidth = 1) +
  geom_line(
    aes(group = Dataset, color = Dataset),
    linewidth = 0.6,
    alpha = 1
  ) +
  geom_point(
    aes(color = Dataset, shape = Dataset),
    size = 2,
    alpha = 1
  ) +
  scale_shape_manual(
    values = c(
      "Schirmer M(2018)" = 2,
      "Feng Q(2015)" = 3,
      "Zeller G(2014)" = 4,
      "Vogtmann E(2016)" = 5,
      "Jie ZY(2017)" = 15,
      "Liu RX(2017)" = 16,
      "Nielsen HB(2014)" = 17,
      "He Q(2017)" = 18,
      "B40-8+VD13" = 0,
      "T4+F1" = 8
    )
  ) +
  stat_pvalue_manual(
    pairwise_res1,
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  scale_y_continuous(
    limits = c(0, 76),
    breaks = seq(0, 76, by = 10)
  )+
  scale_fill_manual(
    values = c(
      "ADAPT" = "#7b95c6",
      "ADMIC.Basic" = "#DB9D5F",
      "ADMIC" = "#c85e62"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Number of differential taxa",
    fill = "Metric",
    color = "Dataset"
  ) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

p2




library(patchwork)
(p2 | p1)+plot_layout(widths = c(5, 5))


