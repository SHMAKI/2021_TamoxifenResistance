library(tidyverse)

# import ggplot theme
source("src/theme.R")
LEGEND_POSITION <- c(0.8, 0.1)

# import data
file <- read.csv("data/growthdata.csv", header=T)
file <- file %>% dplyr::filter(Week>0)

g <-  ggplot(file, aes(x=Week, y=Value, color=Condition, shape = Condition, fill=Condition)) + 
  geom_point(stat="identity", position=position_dodge(width=0.5), size=pointsize, alpha=1) +
  # +/- standard deviation
  stat_summary(fun.data=function(...) mean_se(..., mult=1),#mean_sdl
               geom='errorbar', width=0.1, aes(color=Condition), position=position_dodge(width=0.5)) +
  # points for mean, using hyphens for point shape
  stat_summary(fun.y=mean, aes(colour=Condition), geom='point', shape='-', size=10, position=position_dodge(width=0.5)) +
  # line connecting means
  stat_summary(fun.y=mean,  geom='line', aes(colour=Condition), lty=2, position=position_dodge(width=0.5))
g <- g + ylab("growth rate / week") + xlab("week") + scale_x_continuous(breaks = 1:10)
g <- g + theme_bw(base_size=15) + ng1 + theme(legend.position = LEGEND_POSITION) 
g1 <- g + scale_color_manual(values = alpha(MyColor, 1)) + scale_fill_manual(values = alpha(MyColor, 1))
g2 <- g + scale_colour_manual(values = c("gray", "black")) + scale_fill_manual(values = alpha(c("gray", "black"), 0.7))


#Wletch's test
ttest_list <- c()
for (i in 1:10) {
  w1_c <- file %>% dplyr::filter(Week == i & Condition == "Ctrl")
  w1_i <- file %>% dplyr::filter(Week == i & Condition == "TAM")
  ttest_list[i] <- t.test(w1_c$Value, w1_i$Value)$p.value
}
ttest_list_adj = p.adjust(ttest_list, method = "BH")

x_diff = 0.25
xend = 2.25
ytop = 4.4
ytops = c(ytop+0.1, ytop-0.1)
d.signif <- data.frame(xstart = 1:10 - x_diff,
                       xend =  1:10 + x_diff,
                       #ystart = c(102, 102, 102, 102),
                       ytop = rep(ytops,5),
                       #yend = c(109, 79, 85, 85),
                       label = c(paste0("p = \n", sprintf("%.1e", ttest_list[1])),
                                 paste0(sprintf("%.1e", ttest_list[2])),
                                 paste0(sprintf("%.1e", ttest_list[3])),
                                 paste0(sprintf("%.1e", ttest_list[4])),
                                 paste0(sprintf("%.1e", ttest_list[5])),
                                 paste0(sprintf("%.1e", ttest_list[6])),
                                 paste0(sprintf("%.1e", ttest_list[7])),
                                 paste0(sprintf("%.1e", ttest_list[8])),
                                 paste0(sprintf("%.1e", ttest_list[9])),
                                 paste0(sprintf("%.1e", ttest_list[10]))
                       ),
                       label.x = 1:10,
                       label.y =  rep(ytops+0.1,5)
)

d.signif$label.y[1] = d.signif$ytop[1] + 0.25

gp <- g1 + geom_text(aes(x = label.x, y = label.y, label = label), d.signif, inherit.aes=FALSE, size=3) +
  geom_segment(aes(x = xstart, xend = xend, y = ytop, yend = ytop), d.signif, inherit.aes=FALSE) +
  ylim(c(min(file$Value),5))

# ggsave(file ="Fig/Fig1b.pdf", plot=gp, width=6, height=4)
