# 清空环境变量
rm(list=ls())

# 加载所需的包
library(readxl)
library(ggplot2)
library(ggExtra)
library(vegan)
library(ggthemes)
# 设置工作路径（使用全路径）
setwd("/Users/mayulin/Desktop/data")
# 加载种水平物种丰度表，并设置列名和分隔符
data <- read_excel("tax.xlsx")
data <- as.data.frame(data)
rownames(data) <- data$`Sample`
data <- data[, -1]
data <- t(data)
View(data)

# 加载分组表group
group <- read_excel("group.xlsx")
group<-c(rep("0W",15),rep("1W",15),rep("3W",15))
group_df <- data.frame(group = group)

# 处理空行
data <- data[rowSums(data) > 0, ]

# 处理空行和缺失值
data <- na.omit(data)

#进行PCoA分析
# vegdist函数，计算距离；method，选择距离类型
# 可选euclidean、manhattan、jaccard
bray <- vegdist(data, method = 'bray')
# euclidean <- vegdist(data, method = 'euclidean')
# manhattan <- vegdist(data, method = 'manhattan')
# jaccard <- vegdist(data, method = 'jaccard')
# 计算加权bray-curtis距离
# dune_bray <- vegdist(data, method="bray", binary=F)
bray <- as.matrix(bray)
write.table(bray,"tax.xlsx",sep = "\t") 

# 对非加权距离进行PCoA分析,k等于3,前三个维度。
# 可选k = (nrow(data) - 1
pcoa <- cmdscale(bray, k = 3, eig = T)

# 提取样本点坐标
pcoa_data <- data.frame(pcoa$points)
# 提取列名，便于后续操作。
pcoa_data$Sample_ID <- rownames(pcoa_data)
names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig
eig_percent <- round(pcoa$eig/sum(eig)*100,1)
# 为样本点坐标添加分组信息
# 重新定义group以匹配pcoa_data的行数
# 这里假设你需要根据实际情况调整分组的分配
# 例如，如果分组信息是按顺序分配的，你可以根据pcoa_data的行数重新生成group
new_length <- nrow(pcoa_data)
group <- c(rep("0W", ceiling(new_length / 3)), 
           rep("1W", ceiling(new_length / 3)), 
           rep("3W", new_length - 2 * ceiling(new_length / 3)))
group <- group[1:new_length]  # 确保长度一致

pcoa_result <- cbind(pcoa_data, group)  # PCoA值和分组信息合并

# 进行置换多元（因素）方差分析（PERMANOVA）
# 基于bray-curtis距离进行
dune.div <- adonis2(data ~ group, data = pcoa_result, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis
# 绘制简单的pcoa图
pcoa_result$group <- factor(pcoa_result$group, levels = c("0W", "1W", "3W"))
p = ggplot(pcoa_result, aes(x=PCoA1, y=PCoA2, color=group)) +
  geom_point(aes(color=group,shape=group),size=5)+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       caption = dune_adonis)  +# 也可用title
  scale_colour_manual(values = c("#2f2f2f", "#972b1d", "#0072B2", "#D55E00", "#009E73", "#CC79A7"))+
  theme(
    legend.position.inside = c(0.9, 0.19),  # 使用 legend.position.inside 替代 legend.position
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    axis.text = element_text(color = "black", size = 10)
  )+
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed")+
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")
# 保存
ggsave('bra_pcoa1.pdf', width = 5, height = 5)
p = p + stat_ellipse(data = pcoa_result, 
                     geom = "polygon", 
                     level = 0.9, 
                     linetype = 2, 
                     linewidth = 0.5, 
                     aes(fill = group), 
                     alpha = 0.3, 
                     show.legend = TRUE,
                     type = "norm") +  # 使用普通协方差矩阵
  scale_fill_manual(values = c("#2f2f2f", "#972b1d", "#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  theme(
    # 将图例位置设置到右下角并向上移动
    legend.position = c(0.1, 0.3),  #（水平移动，上下移动）
    legend.background = element_rect(fill = "white", color = "black"),  # 添加图例边框
    legend.box.background = element_rect(fill = "white", color = "black")  # 确保图例框有边框
  )
p

# 添加边缘的密度图
#install.packages("ggExtra") # 如果尚未安装
# 指定图形的宽度、高度和分辨率
pdf(file = "bray_pcoa.pdf", width = 5, height = 5)
## pdf格式
# pdf('marginal.pdf',width = 5,height = 5)
ggMarginal(
  p,
  type =c('density'),
  margins = 'both',
  size = 3.5,
  groupColour = F,
  groupFill = T
) 
dev.off()
