library(tidyverse)

# import ggplot theme
source("src/theme.R")
LEGEND_POSITION <- c(0.75, 0.1)

# import data
file_S <- read.csv("data/Sphase.csv", header=T)
file_G1 <- read.csv("data/G1phase.csv", header=T)

# Fig.1c S phase --------
g <-  ggplot(file_S, aes(x=Week, y=Value, color=Condition, shape = Condition, fill=Condition)) + 
  geom_point(stat="identity", position=position_dodge(width=0.5), size=pointsize, alpha=1) +
  # +/- standard deviation
  stat_summary(fun.data=function(...) mean_se(..., mult=1),#mean_sdl
               geom='errorbar', width=0.1, aes(color=Condition), position=position_dodge(width=0.5)) +
  # points for mean, using hyphens for point shape
  stat_summary(fun=mean, aes(colour=Condition), geom='point', shape='-', size=10, position=position_dodge(width=0.5)) +
  # line connecting means
  stat_summary(fun=mean,  geom='line', aes(colour=Condition), lty=2, position=position_dodge(width=0.5))
g <- g + ylab("S phase (%)") + xlab("week") + scale_x_continuous(breaks = 1:10)
g <- g + theme_bw(base_size=15) + ng1 + theme(legend.position = LEGEND_POSITION) + ylim(c(min(file_S$Value), 45))
g1 <- g + scale_color_manual(values = alpha(MyColor, 1)) + scale_fill_manual(values = alpha(MyColor, 1))
g2 <- g + scale_colour_manual(values = c("gray", "black")) + scale_fill_manual(values = alpha(c("gray", "black"), 0.7))

# Wletch'S test
ttest_list <- c()
for (i in 1:10) {
  w1_c <- file_S %>% dplyr::filter(Week == i & Condition == "Ctrl")
  w1_i <- file_S %>% dplyr::filter(Week == i & Condition == "TAM")
  ttest_list[i] <- t.test(w1_c$Value, w1_i$Value)$p.value
}

ttest_list_adj = p.adjust(ttest_list, method = "BH")

x_diff = 0.25
ytop = 40
ytops = c(ytop-2, ytop+2)
d.signif <- data.frame(xstart = 1:10 - x_diff,
                       xend =  1:10 + x_diff,
                       ytop = rep(ytops,5),
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
                       label.y =  rep(ytops+1,5)
)
d.signif$label.y[1] = d.signif$ytop[1] + 2

gp <- g1 + geom_text(aes(x = label.x, y = label.y, label = label), d.signif, inherit.aes=FALSE, size=3) +
  geom_segment(aes(x = xstart, xend = xend, y = ytop, yend = ytop), d.signif, inherit.aes=FALSE)

# ggsave(file ="Fig/Fig1c.pdf", plot=gp, width=6, height=4)


# Fig.1d G1 phase --------
LEGEND_POSITION <- c(0.8, 0.9)

g <-  ggplot(file_G1, aes(x=Week, y=Value, color=Condition, shape = Condition, fill=Condition)) + 
  geom_point(stat="identity", position=position_dodge(width=0.5), size=pointsize, alpha=1) +
  # +/- standard deviation
  stat_summary(fun.data=function(...) mean_se(..., mult=1),#mean_sdl
               geom='errorbar', width=0.1, aes(color=Condition), position=position_dodge(width=0.5)) +
  # points for mean, using hyphens for point shape
  stat_summary(fun=mean, aes(colour=Condition), geom='point', shape='-', size=10, position=position_dodge(width=0.5)) +
  # line connecting means
  stat_summary(fun=mean,  geom='line', aes(colour=Condition), lty=2, position=position_dodge(width=0.5))
g <- g + ylab("G1 phase (%)") + xlab("week") + scale_x_continuous(breaks = 1:10)
g <- g + theme_bw(base_size=15) + ng1 + theme(legend.position = LEGEND_POSITION) + ylim(min(file_G1$Value),90) 
#+ ylim(c(min(ds$mean-ds$sd),1.8)) #+ scale_shape_manual(values=c(16,17))

g1 <- g + scale_color_manual(values = alpha(MyColor, 1)) + scale_fill_manual(values = alpha(MyColor, 1))
g2 <- g + scale_colour_manual(values = c("gray", "black")) + scale_fill_manual(values = alpha(c("gray", "black"), 0.7))

# Wletch'S test
ttest_list2 <- c()
for (i in 1:10) {
  w1_c <- file_G1 %>% dplyr::filter(Week == i & Condition == "Ctrl")
  w1_i <- file_G1 %>% dplyr::filter(Week == i & Condition == "TAM")
  ttest_list2[i] <- t.test(w1_c$Value, w1_i$Value)$p.value
}
ttest_list_adj = p.adjust(ttest_list2, method = "BH")

x_diff = 0.25
ytop = 85
ytops = c(ytop+2, ytop-2)
d.signif2 <- data.frame(xstart = 1:10 - x_diff,
                       xend =  1:10 + x_diff,
                       ytop = rep(ytops,5),
                       label = c(paste0("p = \n", sprintf("%.1e", ttest_list2[1])),
                                 paste0(sprintf("%.1e", ttest_list2[2])),
                                 paste0(sprintf("%.1e", ttest_list2[3])),
                                 paste0(sprintf("%.1e", ttest_list2[4])),
                                 paste0(sprintf("%.1e", ttest_list2[5])),
                                 paste0(sprintf("%.1e", ttest_list2[6])),
                                 paste0(sprintf("%.1e", ttest_list2[7])),
                                 paste0(sprintf("%.1e", ttest_list2[8])),
                                 paste0(sprintf("%.1e", ttest_list2[9])),
                                 paste0(sprintf("%.1e", ttest_list2[10]))
                       ),
                       label.x = 1:10,
                       label.y =  rep(ytops+1,5)
)

d.signif2$ytop[1] = d.signif2$ytop[1] - 2
d.signif2$label.y[1] = d.signif2$ytop[1] + 2
d.signif2$ytop[4:10] = d.signif2$ytop[4:10] - 10
d.signif2$label.y[4:10] = d.signif2$label.y[4:10] - 10

gp <- g1 + geom_text(aes(x = label.x, y = label.y, label = label), d.signif2, inherit.aes=FALSE, size=3) +
  geom_segment(aes(x = xstart, xend = xend, y = ytop, yend = ytop), d.signif2, inherit.aes=FALSE)

# ggsave(file ="Fig/Fig1d.pdf", plot=gp, width=6, height=4)
