rm(list = ls())
# 加载所需的包
library(agricolae)
library(ggplot2)
library(grid)
library(vegan)
library(readxl)
library(dplyr)

# 设置工作路径
setwd("/Users/xinhao/Desktop/test")

# 读取物种丰度数据
tax <- read.csv("chao1.boxplot_0.csv")
# 读取分组信息
group <- read_excel("group.xlsx")

tax <- as.data.frame(tax)
tax <- tax[,-1]
tax <- t(tax)
tax <- as.data.frame(tax)

tax$sample <- rownames(tax)
colnames(group)[1] <- "sample"
colnames(group)[2] <- "group"
tax <- merge(tax, group, by = "sample")

# 明确选择物种丰度列（排除'sample'和'group'）
species_columns <- setdiff(names(tax), c("sample", "group"))
species_data <- tax[, species_columns]
species_data <- as.data.frame(lapply(species_data, as.numeric))

# 处理非整数类型的数据，将其转换为整数
species_data <- as.data.frame(lapply(species_data, round))

alpha <- alpha_index(tax_numeric)
alpha$name <- rownames(alpha) 
alpha$group <- c(rep('0W', 15), rep('1W', 15),rep('3W', 15))
write.csv(alpha, "alpha_index.csv", row.names = T )#将所有alpha指数写入一个文档中

# 计算指数
#tax$shannon_entropy <- diversity(species_data,index = "shannon")
#tax$shannon_entropy <- diversity(species_data, index = "simpson")
#tax$shannon_entropy <- estimateR(species_data)["S.chao1", ]
tax$shannon_entropy <- estimateR(species_data)["S.ACE", ]

# 执行ANOVA分析
anova_result <- aov(simpson_entropy ~ group, data = tax)
# 执行Tukey HSD检验
tukey_result <- TukeyHSD(anova_result)
index <- as.data.frame(tukey_result$groups)

# 重新设置group因子的水平顺序
tax$group <- factor(tax$group, levels = c("0W", "1W", "3W"))

# 计算每个组的最大值并添加Tukey检验标签
s2 <- tax %>%
  group_by(group) %>%
  summarise(shannon_entropy = max(shannon_entropy))
s2$label <- index$groups[match(rownames(index), s2$group)]

# 绘制alpha多样性图

# 自定义 6 种辨识度高的颜色
custom_colors <- c("#6B7280", "#F97316", "#0000FF", "#F97316", "#FF00FF", "#C1A57B")

p1 <- ggplot(tax, aes(x = group, y = value, fill = group, color = group)) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 1) +
  geom_boxplot(outlier.shape = NA,alpha = 0.6, width = 0.8, size = 1.5) +
  # 使用手动指定的填充颜色
  scale_fill_manual(values = custom_colors) +
  # 使用手动指定的边框颜色
  scale_color_manual(values = custom_colors) +
  labs(y = "Chao1 richness", x = "") +
  geom_signif(comparisons = list(c("0W", "1W"), c("0W", "3W"), c("1W", "3W")),
             annotations = c("ns", "***", "ns"),  # 根据p值设置：* p<0.05, ** p<0.01, *** p<0.001
             y_position = c(120, 160, 135),
             size = 0.8,              # 较粗的标线
             textsize = 5,          # 较大的文本
             fontface = "bold",     # 粗体文本
             family = "serif" ,
             color = 'black',
             tip_length = 0.03)+
  theme(
    panel.grid = element_line(color = 'white', size = 2),
    panel.background = element_rect(size = 1, fill = '#EBEBEB', color = 'transparent'),
    panel.border = element_rect(fill = 'transparent', color = 'transparent'),
    axis.title.y = element_text(family = "serif", size = 20),
    axis.line.y = element_blank(),
    axis.text.y = element_text(family = "serif", size = 15, color = 'black'),
    axis.text.x = element_text(family = "serif", size = 15, color = 'black'),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.ticks.length.y = unit(0.25, 'cm'),
    strip.text = element_text(face = 'bold', family = "serif", size = 25),
    strip.background = element_rect(fill = 'gray85', color = 'transparent'),
    legend.key.size = unit(1, "cm"),           # 适中的图例键大小
    legend.text = element_text(size = 14),     # 稍大的图例文本
    legend.title = element_text(size = 16))  # 突出显示的图例标题                            )
  
p1
# 保存图形为PDF格式
ggsave("Chao1.pdf", plot = p1, width = 8, height = 6, dpi = 300)



