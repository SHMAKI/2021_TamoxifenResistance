library(ggplot2)
library(Hmisc)
rm(list=ls()) #全オブジェクトの削除

file <- read.csv("growthdata.csv", header=T)

p <- ggplot(aes(x=Week, y=Value, color=Condition), data=file)
p <- p + geom_point() #+ ylim(0, 1.1)
plot(p)

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, geom=geom, linewidth=2, ...)
}
stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, geom=geom, size = 0.5,  color=factor(Condition), ...)
}

q <- qplot(Week, Value, data=file)  + xlab("Week") + ylab("Growth rate / week") + ylim(0, 4.5) + scale_x_continuous(breaks=0:10)+ theme_bw(base_size = 18) + theme(legend.position = "bottom", legend.direction = "horizontal")
q1 <- q + stat_summary(fun.y = mean, geom = "line", aes(color=factor(Condition))) + geom_point(aes(color=factor(Condition)))
#q2 <- q + stat_sum_single(mean, geom="line", aes(color=factor(Condition))) 
#q3 <- q + stat_sum_single(mean, geom="line", colour="black") + aes(color=factor(Run))

ggsave(file ="./1.pdf", plot=q1, dpi = 300, width=16, height=9)
ggsave(file ="./2.pdf", plot=q2, dpi = 300, width=16, height=9)
ggsave(file ="./3.pdf", plot=q3, dpi = 300, width=16, height=9)
ggsave(file ="./4.pdf", plot=qq2, dpi = 300, width=16, height=9)