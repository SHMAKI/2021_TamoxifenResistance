library(tidyverse)
library(plyr)

# import ggplot theme
source("src/theme.R")

# chemical_condition
MyLabels <- c("Ctrl",expression("GSK467 10 "~mu~"M"))

# experimental_condition
LABEL1 <- c("W0","W3","W6","W9")

# import data
file <- read_csv("data/MTTassay.csv")
file <- file %>% mutate(Condition = paste0("W",week))
file$chemical[file$chemical == "GSK467_10"] <- "GSK467 10 µM"
file$Condition <- factor(file$Condition, levels=LABEL1)
file$chemical <- factor(file$chemical, levels=c("Ctrl","GSK467 10 µM"))
file <- with(file, file[order(Condition),])
ds <- ddply(file, .(Condition, si, chemical), summarise, MEAN = mean(value), SD = sd(value))

# Wletch's test at W0, W3, W6, W9
w0_c <- file %>% dplyr::filter(week == 0 & si == "siCtrl" & chemical == "Ctrl")
w0_i <- file %>% dplyr::filter(week == 0 & si == "siPML" & chemical == "GSK467 10 µM")
w0_onlyPML <- file %>% dplyr::filter(week == 0 & si == "siPML" & chemical == "Ctrl")
w0_onlyGSK <- file %>% dplyr::filter(week == 0 & si == "siCtrl" & chemical == "GSK467 10 µM")
w0_wel <- t.test(w0_c$value, w0_i$value)
w0_wel_PMLvsBoth <- t.test(w0_onlyPML$value, w0_i$value)
w0_wel_GSKvsBoth <- t.test(w0_onlyGSK$value, w0_i$value)

w3_c <- file %>% dplyr::filter(week == 3 & si == "siCtrl" & chemical == "Ctrl")
w3_i <- file %>% dplyr::filter(week == 3 & si == "siPML" & chemical == "GSK467 10 µM")
w3_onlyPML <- file %>% dplyr::filter(week == 3 & si == "siPML" & chemical == "Ctrl")
w3_onlyGSK <- file %>% dplyr::filter(week == 3 & si == "siCtrl" & chemical == "GSK467 10 µM")
w3_wel <- t.test(w3_c$value, w3_i$value)
w3_wel_PMLvsBoth <- t.test(w3_onlyPML$value, w3_i$value)
w3_wel_GSKvsBoth <- t.test(w3_onlyGSK$value, w3_i$value)

w6_c <- file %>% dplyr::filter(week == 6 & si == "siCtrl" & chemical == "Ctrl")
w6_i <- file %>% dplyr::filter(week == 6 & si == "siPML" & chemical == "GSK467 10 µM")
w6_onlyPML <- file %>% dplyr::filter(week == 6 & si == "siPML" & chemical == "Ctrl")
w6_onlyGSK <- file %>% dplyr::filter(week == 6 & si == "siCtrl" & chemical == "GSK467 10 µM")
w6_wel <- t.test(w6_c$value, w6_i$value)
w6_wel_PMLvsBoth <- t.test(w6_onlyPML$value, w6_i$value)
w6_wel_GSKvsBoth <- t.test(w6_onlyGSK$value, w6_i$value)

w9_c <- file %>% dplyr::filter(week == 9 & si == "siCtrl" & chemical == "Ctrl")
w9_i <- file %>% dplyr::filter(week == 9 & si == "siPML" & chemical == "GSK467 10 µM")
w9_onlyPML <- file %>% dplyr::filter(week == 9 & si == "siPML" & chemical == "Ctrl")
w9_onlyGSK <- file %>% dplyr::filter(week == 9 & si == "siCtrl" & chemical == "GSK467 10 µM")
w9_wel <- t.test(w9_c$value, w9_i$value)
w9_wel_PMLvsBoth <- t.test(w9_onlyPML$value, w9_i$value)
w9_wel_GSKvsBoth <- t.test(w9_onlyGSK$value, w9_i$value)

# adj multiple test
W0_FDR <- p.adjust(c(w0_wel$p.value, w0_wel_PMLvsBoth$p.value, w0_wel_GSKvsBoth$p.value), "BH")
W3_FDR <- p.adjust(c(w3_wel$p.value, w3_wel_PMLvsBoth$p.value, w3_wel_GSKvsBoth$p.value), "BH")
W6_FDR <- p.adjust(c(w6_wel$p.value, w6_wel_PMLvsBoth$p.value, w6_wel_GSKvsBoth$p.value), "BH")
W9_FDR <- p.adjust(c(w9_wel$p.value, w9_wel_PMLvsBoth$p.value, w9_wel_GSKvsBoth$p.value), "BH")

# coordinates for ctrl vs both
x <- 0.75
xend <- 2.25 + 0.1
ytop <- c(125+15, 105+10, 105+10, 105+10)
d.signif_revise1 <- data.frame(x = c(x,x,x,x),
                       xend = c(xend, xend, xend, xend),
                       ystart = c(102, 102, 102, 102),
                       ytop = ytop,
                       yend = c(109, 79, 85, 85),
                       label = c(paste0("q = ", sprintf("%.1e", W0_FDR[1])),
                                 paste0(sprintf("%.1e", W3_FDR[1])),
                                 paste0(sprintf("%.1e", W6_FDR[1])),
                                 paste0(sprintf("%.1e", W9_FDR[1]))
                       ),
                       label.x = c(1.5, 1.5, 1.5, 1.5),
                       label.y = ytop+c(3,2,2,2),
                       Condition = c("W0", "W3", "W6", "W9")
)

# coordinates for GSKonly vs both
x = 1.25
xend = 2.25
ytop = c(125+10, 105+5, 105+5, 105+5)
d.signif_revise2 <- data.frame(x = c(x,x,x,x),
                               xend = c(xend, xend, xend, xend),
                               ystart = c(102, 85, 88, 90),
                               ytop = ytop,
                               yend = c(109, 79, 85, 85),
                               label = c(paste0("q = ", sprintf("%.1e", W0_FDR[3])),
                                         paste0(sprintf("%.1e", W3_FDR[3])),
                                         paste0(sprintf("%.1e", W6_FDR[3])),
                                         paste0(sprintf("%.1e", W9_FDR[3]))
                               ),
                               label.x = c(1.75, 1.75, 1.75, 1.75),
                               label.y = ytop+c(3,2,2,2),
                               Condition = c("W0", "W3", "W6", "W9")
)

# coordinates for siPMLonly vs both
x = 1.75
xend = 2.25 - 0.1
ytop = c(125, 105, 105, 105)
d.signif_revise3 <- data.frame(x = c(x,x,x,x),
                               xend = c(xend, xend, xend, xend),
                               ystart = c(120, 92, 92, 92),
                               ytop = ytop,
                               yend = c(109, 79, 85, 85),
                               label = c(paste0("q = \n", sprintf("%.1e", W0_FDR[2])),
                                         paste0(sprintf("%.1e", W3_FDR[2])),
                                         paste0(sprintf("%.1e", W6_FDR[2])),
                                         paste0(sprintf("%.1e", W9_FDR[2]))
                               ),
                               label.x = c(2, 2, 2, 2)-0.2,
                               label.y = ytop+c(5,2,2,2),
                               Condition = c("W0", "W3", "W6", "W9")
)


# draw graph
g <-  ggplot(file, aes(x=si, y=value, color=chemical, shape = chemical, fill=chemical)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color="black", width=0.5, size=LINEWIDTH/4, position=position_dodge(width=0.99), lty=1) +
  geom_point(stat="identity", position=position_dodge(width=0.99), size=pointsize, alpha=0.5) + 
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y.., color=chemical), geom='errorbar', size=LINEWIDTH/2, position=position_dodge(width=0.99))

g <- g + theme_bw(base_size=15) + ng1 + theme(axis.title.x=element_blank())
g_2 <- g + facet_grid(.~Condition) + ylab("Relative absorbance (%)")
g1 <- g_2 + scale_colour_manual(values = MyColor) + scale_fill_manual(values = alpha(MyColor, 1))

gp_revise <- g1 +  geom_text(aes(x = label.x, y = label.y, label = label), d.signif_revise1, inherit.aes=FALSE, ) +
  geom_segment(aes(x = x, xend = x, y = ystart , yend = ytop), d.signif_revise1, inherit.aes=FALSE) +
  geom_segment(aes(x = x, xend = xend, y = ytop, yend = ytop), d.signif_revise1, inherit.aes=FALSE) +
  geom_segment(aes(x = xend, xend = xend, y = ytop, yend = yend), d.signif_revise1, inherit.aes=FALSE)
gp_revise <- gp_revise +  geom_text(aes(x = label.x, y = label.y, label = label), d.signif_revise2, inherit.aes=FALSE, ) +
  geom_segment(aes(x = x, xend = x, y = ystart , yend = ytop), d.signif_revise2, inherit.aes=FALSE) +
  geom_segment(aes(x = x, xend = xend, y = ytop, yend = ytop), d.signif_revise2, inherit.aes=FALSE) +
  geom_segment(aes(x = xend, xend = xend, y = ytop, yend = yend), d.signif_revise2, inherit.aes=FALSE)
gp_revise <- gp_revise +  geom_text(aes(x = label.x, y = label.y, label = label), d.signif_revise3, inherit.aes=FALSE, ) +
  geom_segment(aes(x = x, xend = x, y = ystart , yend = ytop), d.signif_revise3, inherit.aes=FALSE) +
  geom_segment(aes(x = x, xend = xend, y = ytop, yend = ytop), d.signif_revise3, inherit.aes=FALSE) +
  geom_segment(aes(x = xend, xend = xend, y = ytop, yend = yend), d.signif_revise3, inherit.aes=FALSE)

ggsave(file ="Fig/fig5a.pdf", plot=gp_revise, width=8, height=5)