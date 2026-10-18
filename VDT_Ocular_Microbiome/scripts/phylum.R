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
data <- read_excel("phylum.xlsx")
View(data)
# 找出重复的行
duplicated(data$phylum)
duplicated_rows <- duplicated(data$phylum)

# 去除重复的行
data <- data[!duplicated_rows, ]
data<-as.data.frame(data)

# 检查 phylum 列是否有缺失值并替换
if (any(is.na(data$phylum))) {
  data$phylum[is.na(data$phylum)] <- "MissingPhylum"  # 用"MissingPhylum"替换缺失值
}

# 再尝试将 phylum 列设置为行名
rownames(data)<-data$phylum
# 提取 taxonomy 列
taxonomy_column <- data$tax

# 删除第一列（根据实际需求确认是否需要删除）
data <- data[,-1]

# 检查每一列是否可以转换为数值型
can_convert <- sapply(data, function(x) all(is.na(as.numeric(x)) | is.numeric(x)))
# 选择可以转换的列进行转换
data[, can_convert] <- apply(data[, can_convert], 2, as.numeric)

# 根据行总和降序排序并取前 10 行
top_10_indices <- order(rowSums(data), decreasing = TRUE)[1:10]
taxonomy_column <- taxonomy_column[top_10_indices]
data <- data[top_10_indices, ]

# 查看数据
head(data)

# 检查数据框中的缺失值
if (any(is.na(data))) {
  # 用 0 填充缺失值（你可以根据实际情况选择其他填充值）
  data[is.na(data)] <- 0
}

# 删除 tax 列
#data <- data[,!(colnames(data) == "tax")]

# 删除列和为 0 的列
data <- data[, colSums(data) != 0]
print(colSums(data))
# 查看数据
print(data)

# 对每列进行归一化处理
for(i in 1:ncol(data)){
  data[,i] <- data[,i]/colSums(data)[i]
}

# View(data_normalized)
# data<-data_normalized
colnames(data)
library(reshape2)
data$tax<-tax
data_long <- melt(data)
View(data_long)
data_long$group <- sub("_.*", "", data_long$variable)

# 查看结果
head(data_long)
View(data_long)
data_long$group<-substr(data_long$group,1,2)
p2<- ggplot(data = data_long, aes(x = variable, y = value, fill = tax,
                                  stratum = tax, alluvium = tax)) +
  geom_col(position ='stack', width = 0.6)+  # stack：堆叠图
  geom_stratum(width = 0.6, color='white')+# 添加连接线
  geom_alluvium(alpha = 0.4,
                width = 0.6,
                color='white',
                linewidth = 1,
                curve_type = "linear")+
  scale_fill_manual(
    values =  rev(c(
      "#2f2f2f", "#972b1d", "#0072B2", "#00A087B2",
      "#F39B7FB2", "#91D1C2B2", "#8491B4B2", "#DC0000B2",
      "#7E6148B2", "yellow", "darkolivegreen1", "lightskyblue",
      "orange", "purple", "brown", "pink", "gray", "cyan",
      "magenta", "chartreuse", "deeppink", "#999999"
    ))
  )+
  labs(x = 'Samples', y = 'Relative Abundance(%)') +
  facet_wrap(~group,scales = 'free_x', ncol = 3)+
  theme(strip.text = element_text(color ="black", size = 12),# 自定义分面文本
        strip.background = element_rect(color = "black", fill="grey90"))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        strip.text = element_text(size = 12)) + #设置主题背景，根据自己的需求定制
  # angle：调整横轴标签倾斜角度
  # hjust：上下移动横轴标签
  theme(axis.text.x=element_text(angle=45, hjust=1),axis.text = element_text(size = 12),
        axis.title = element_text(size = 13), legend.title = element_blank(),
        legend.text = element_text(size = 11),legend.position="bottom",
        legend.key = element_rect(fill = NA, color = NA))+guides(fill=guide_legend(ncol = 5,byrow = T,
                                                                                   override.aes = list(size = 3),title.position = "top"))
p2
ggsave("phylum.pdf",p2,width=12,height=6,family="serif")


