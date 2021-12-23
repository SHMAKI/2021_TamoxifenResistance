#cor_between cell-to-cell in each time point
library(amap)
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(grid)

#from monocle dat, cds
unique(Type)
SD_rank <- apply(umi_matrix,1,sd) #/apply(dat5,1,mean)
TOP <- nrow(umi_matrix) + 1 #20000 #5001 1001
dat_SD_top_1000 <- umi_matrix [rank(-SD_rank, na.last = TRUE) < TOP,]

dat5_CTRL <- dat_SD_top_1000[,grep("Control", colnames(dat_SD_top_1000))]
dat5_W3 <- dat_SD_top_1000[,grep("TAM3weeks", colnames(dat_SD_top_1000))]
dat5_W6 <- dat_SD_top_1000[,grep("TAM6weeks", colnames(dat_SD_top_1000))]
dat5_W9 <- dat_SD_top_1000[,grep("TAM9weeks", colnames(dat_SD_top_1000))]

dat5_CTRL_cor <- cor(dat5_CTRL, method = "p")
dat5_W3_cor <- cor(dat5_W3, method = "p")
dat5_W6_cor <- cor(dat5_W6, method = "p")
dat5_W9_cor <- cor(dat5_W9, method = "p")
#X[upper.tri(X)]で上三角成分
dat5_CTRL_cor_v <- dat5_CTRL_cor[upper.tri(dat5_CTRL_cor)]
dat5_W3_cor_v <- dat5_W3_cor[upper.tri(dat5_W3_cor)]
dat5_W6_cor_v <- dat5_W6_cor[upper.tri(dat5_W6_cor)]
dat5_W9_cor_v <- dat5_W9_cor[upper.tri(dat5_W9_cor)]

DAT_cor <- data.frame(value=c(dat5_CTRL_cor_v, dat5_W3_cor_v, dat5_W6_cor_v, dat5_W9_cor_v),
                      week=c(rep("0", length(dat5_CTRL_cor_v)), 
                              rep("3", length(dat5_W3_cor_v)),
                              rep("6", length(dat5_W6_cor_v)),
                              rep("9", length(dat5_W9_cor_v))
                              )
                      )
#自己相関以外の相関係数の平均値
dat5_CTRL_cor_v2 <- apply(dat5_CTRL_cor, 1, function(x) {(sum(x) - 1) / (length(x) -1)})
dat5_W3_cor_v2 <-  apply(dat5_W3_cor, 1, function(x) {(sum(x) - 1) / (length(x) -1)})
dat5_W6_cor_v2 <-  apply(dat5_W6_cor, 1, function(x) {(sum(x) - 1) / (length(x) -1)})
dat5_W9_cor_v2 <-  apply(dat5_W9_cor, 1, function(x) {(sum(x) - 1) / (length(x) -1)})

DAT_cor_mean <- data.frame(value=c(dat5_CTRL_cor_v2, dat5_W3_cor_v2, dat5_W6_cor_v2, dat5_W9_cor_v2),
                           week=c(rep("0", length(dat5_CTRL_cor_v2)), 
                              rep("3", length(dat5_W3_cor_v2)),
                              rep("6", length(dat5_W6_cor_v2)),
                              rep("9", length(dat5_W9_cor_v2))
                      )
)




textsize <- 18
pointsize <- 2.5
legendsize <- 18
LINEWIDTH <- 2
LEGEND_POSITION <- c(0.8, 0.9)
ng1 <- theme(#panel.background = element_rect(fill = "white",colour = "white"),
  #panel.grid.major = element_line(colour = NA),
  #axis.line = element_line(size = 1.2, colour="black"),
  #axis.ticks.x=element_blank(),
  #             axis.ticks=element_line(color="black"),
  #panel.grid.minor = element_line(colour = NA),
  axis.text=element_text(color="black",size=textsize),
  axis.title=element_text(color="black",size=textsize),
  #axis.line = element_line(size=2),
  panel.border = element_rect(color = "black", fill = NA, size = LINEWIDTH),
  panel.background = element_rect(color = "black", fill = NA, size = LINEWIDTH),
  
  #panel.margin=unit(.05, "lines"),
  #legend.position = "bottom", #LEGEND_POSITION,#"right",# c(0.7, 1),   
  #legend.direction = "horizontal", #or "horizontal",
  legend.background=element_rect(fill="transparent"),
  legend.key=element_rect(linetype = 0, fill = "white"),
  legend.text = element_text(size=legendsize),
  legend.title = element_text(size=legendsize), #element_blank(),
  legend.key.size = unit(1, "cm"),
  strip.background = element_rect(color="transparent", fill="transparent", size=LINEWIDTH),
  strip.text = element_text(size=textsize)             
)

DIR8 <- "results8"
dir.create(DIR8)
g <- ggplot(DAT_cor, aes(x=week, y=value)) + geom_violin(fill="white") +
  geom_boxplot(width=.05, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") + 
  #stat_summary(fun.y=median, geom="point", fill="white", shape=21) +
  theme_bw() + scale_fill_aaas() + ng1 + ylab("Correlation coefficient\nbetween cells") + xlab("week") + ylim(c(0.8, 1.0))#scale_y_log10() +
  ggsave(file =paste0(DIR8, "/Fig_cor_in_each_timepoint_top", TOP-1, "genes.pdf"), plot=g, width=4, height=4)

g <- ggplot(DAT_cor, aes(x=week, y=value)) + #geom_violin(fill="white") +
    geom_boxplot(width=.5, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") + #scale_y_log10() +
    #stat_summary(fun.y=median, geom="point", fill="white", shape=21) +
    theme_bw() + ng1 + ylab("Correlation coefficient\nbetween cells") + xlab("week") #scale_y_log10() +
my_comparisons <- list( c("0", "3"), c("3", "6"), c("6", "9") ) # Add pairwise comparisons p-value
g <- g + stat_compare_means(comparisons = my_comparisons, bracket.size = 0.5, step.increase = 0.02, tip.length = 0.01) #+ stat_compare_means(label.y = 15)     # Add global p-value
#g <- g + coord_cartesian(ylim=c(0.8, 1.1))
#g <- g + ylim(c(0.8, 1.1))
#ggsave(file =paste0(DIR8, "/Fig_cor_in_each_timepoint_top_ylim.pdf"), plot=g, width=4, height=4)
#ggsave(file =paste0(DIR8, "/Fig_cor_in_each_timepoint_top", TOP-1, "gene_revised_boxplot.pdf"), plot=g, width=4, height=4)
  

g2 <- ggplot(DAT_cor_mean, aes(x=week, y=value)) + geom_violin(fill="white") +
    geom_boxplot(width=.05, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") + 
    #stat_summary(fun.y=median, geom="point", fill="white", shape=21) +
    theme_bw() + scale_fill_aaas() + ng1 + ylab("Correlation coefficient between cells") + xlab("week") + ylim(c(0.8, 1.0))#scale_y_log10() +
  ggsave(file =paste0(DIR8, "/Fig_cor_in_each_timepoint_top", TOP-1, "genes_mean.pdf"), plot=g2, width=4, height=4)
