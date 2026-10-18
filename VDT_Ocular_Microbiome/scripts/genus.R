# 清除工作空间中的所有对象
rm(list = ls())
# 设置工作路径
setwd("/Users/mayulin/Desktop/data")
# 加载所需的包
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggalluvial)
library(ggh4x)
library(tibble)
library(dplyr)
library(readxl)

# 加载属水平物种丰度表（确保行名为菌种名称，列为样本）
# 读取数据
# 读取数据
data <- read_excel("genus.xlsx")

# 找出重复的行
duplicated_rows <- duplicated(data[[1]])
# 去除重复的行
data <- data[!duplicated_rows, ]

# 检查第一列是否有缺失值并替换
if (any(is.na(data[[1]]))) {
  data[[1]][is.na(data[[1]])] <- "MissingPhylum"  # 用"MissingPhylum"替换缺失值
}

# 移除现有的行名
rownames(data) <- NULL

# 再尝试将第一列设置为行名
data <- data %>%
  column_to_rownames(var = colnames(.)[1])

# 计算每一行的和并按降序排列
data$sum <- rowSums(data, na.rm = TRUE)
data <- data[order(-data$sum), ]
data$sum <- NULL  # 移除sum列

# 过滤掉Unclassified行
#data <- data[!grepl("s_Unclassified", rownames(data)), ]

# 计算分组平均值（根据实际样本分组调整列索引）
data_avg <- data.frame(
  "0W" = rowMeans(data[, 1:16]),    # 调整列索引为NL组的样本
  "1W" = rowMeans(data[, 17:32]), 
  "3W" = rowMeans(data[, 33:48])
)

# 按总丰度排序并取前20
total <- colSums(data_avg)
data_top <- data_avg[order(rowSums(data_avg), decreasing = TRUE)[1:20], ]

# 计算相对丰度并添加Others
data_rel <- data_top / total
others <- 1 - colSums(data_rel)
data_rel <- rbind(data_rel, Others = others)

# 准备绘图数据
plot_data <- data_rel %>%
  rownames_to_column("Taxonomy") %>%
  melt(id.vars = "Taxonomy") %>%
  mutate(Taxonomy = factor(Taxonomy, levels = rev(c(rownames(data_rel)[1:20], "Others")))) %>%
  filter(!Taxonomy %in% c("Unclassified", "Others"))

# 过滤掉 Taxonomy 列中包含 "NA" 的行
plot_data_filtered <- plot_data %>% 
  filter(!grepl("NA", Taxonomy))

# 绘图
p <- ggplot(plot_data_filtered, aes(x = variable, y = 100 * value, fill = Taxonomy)) +
  geom_col(position = "stack", width = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "#2f2f2f", "#972b1d", "#0072B2", "#00A087B2", 
    "#F39B7FB2", "#91D1C2B2", "#8491B4B2", "#DC0000B2", 
    "#7E6148B2", "yellow", "darkolivegreen1", "lightskyblue",
    "orange", "purple", "brown", "pink", "gray", "cyan", 
    "magenta", "chartreuse", "deeppink", "#999999"
  )) +
  labs(x = "", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 1))

# 显示图形
p
# 保存图形
ggsave("genus.pdf", plot = p, width = 12, height = 9, dpi = 600)

# 过滤掉 Taxonomy 列中包含 "NA" 的行
plot_data_filtered <- plot_data %>% 
  filter(!grepl("NA", Taxonomy))

# 绘图
p <- ggplot(plot_data_filtered, aes(x = variable, y = 100 * value, fill = Taxonomy,
                                    stratum = Taxonomy, alluvium = Taxonomy)) +
  geom_col(position ='stack', width = 0.6) +# width控制柱子宽度，设置为1可紧挨着
  geom_stratum(width = 0.6, color='white')+# 添加连接线
  geom_alluvium(alpha = 0.4,
                width = 0.6,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values =  rev(c(
    "#2f2f2f", "#972b1d", "#0072B2",  "#00A087B2", 
    "#F39B7FB2", "#91D1C2B2", 
    "#8491B4B2", "#DC0000B2", 
    "#7E6148B2", "yellow", 
    "darkolivegreen1", "lightskyblue",
    "orange", "purple", "brown", 
    "pink", "gray", "cyan", "magenta",
    "chartreuse", "deeppink"))) +
  labs(x = '', y = 'Relative Abundance (%)') +
  # theme_minimal() +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12),
        axis.title = element_text(size = 13), legend.title = element_blank(), 
        legend.text = element_text(size = 11), legend.position = "right",
        legend.key = element_rect(fill = NA, color = NA)) + guides(fill=guide_legend(ncol=1,byrow=TRUE))

# 显示图形
p
##保存图形
ggsave("genus1.pdf", plot = p, width = 12, height = 9, dpi = 600)

