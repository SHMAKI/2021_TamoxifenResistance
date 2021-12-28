textsize <- 18
pointsize <- 2.5
legendsize <- 12
LINEWIDTH <- 2
MyColor <- c("#CC6677", "skyblue")

ng1 <- theme(
  axis.text=element_text(color="black",size=textsize),
  axis.title=element_text(color="black",size=textsize),
  panel.border = element_rect(color = "black", fill = NA, size = LINEWIDTH),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.background=element_rect(fill="transparent"),
  legend.key=element_rect(linetype = 0, fill = "white"),
  legend.text = element_text(size=legendsize),
  legend.title = element_blank(),legend.key.size = unit(1, "cm"),
  strip.background = element_rect(color="transparent", fill="transparent", size=LINEWIDTH),
  strip.text = element_text(size=textsize)             
)