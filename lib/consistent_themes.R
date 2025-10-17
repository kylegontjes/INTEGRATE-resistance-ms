# Consistent themes

## Packages
library(gridExtra)
library(ggtext)
library(kableExtra)
library(knitr)

# Manuscript Background theme
theme_bw_me <- theme(panel.background = element_rect(fill = "white",colour = NA), panel.grid = element_blank(),
                     strip.background = element_rect(fill = "white",colour = "black",linewidth=1), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),legend.position = "bottom")

# Study theme
study_color_fill <- scale_fill_manual(breaks = c("2014-15","2021-23"), values = c("blue","red"),name="Study")

# Sequence type
ST_scale_fill <- scale_fill_manual(breaks = c("ST258","Not ST258"),values = c("#D6AE7F",'black'),name = "ST")