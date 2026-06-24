# Consistent themes

## Packages
library(gridExtra)
library(ggtext)
library(kableExtra)
library(knitr)
library(tidyverse)

# Common functions
source("./lib/common_functions.R")

# Manuscript Background theme
theme_bw_me <- theme(panel.background = element_rect(fill = "white",colour = NA), panel.grid = element_blank(),
                     strip.background = element_rect(fill = "white",colour = "black",linewidth=1), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),legend.position = "bottom")

# Favorite kable
favorite_kable <- function (x){
  x %>% kable(., format = "html", table.attr = "style='width:100%;'",
              row.names = F) %>% kable_styling(bootstrap_options = c("striped",
                                                                     "hover", "condensed", "responsive"))
}

# Genotype scale (for phylogeny)
feature_colors <- c(`1` = "black",`0`="white")
feature_scale <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Genotype", guide = guide_legend(nrow=3,order=3, title.position = "top", label.position = "right"))
resistance_bin_scale <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Resistant","Susceptible"),name="Resistance", guide = guide_legend(nrow=3,order=3, title.position = "top", label.position = "right"))
carbapenemase_presence_scale <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Carbapenemase presence", guide = guide_legend(nrow=4,order=3, title.position = "top", label.position = "right"))

# ST scale
ST_palette <- readRDS("./lib/ST_palette.RDS")
ST_breaks <- names(ST_palette)
break_of_sts_renamed <- gsub("Klebsiella pneumoniae ","",ST_breaks)
ST_scale <- scale_fill_manual(breaks =ST_breaks,values = ST_palette,labels = break_of_sts_renamed, name = "Sequence type (ST)", guide = guide_legend(order=1,ncol=7, title.position = "top", label.position = "right"))
ST_scale_3_row <- scale_fill_manual(breaks = ST_breaks,values = ST_palette,labels = break_of_sts_renamed, name = "Sequence type (ST)", guide = guide_legend(order=1,nrow=3, title.position = "top", label.position = "right"))
ST_scale_4_row <- scale_fill_manual(breaks = ST_breaks,values = ST_palette,labels = break_of_sts_renamed, name = "Sequence type (ST)", guide = guide_legend(order=1,nrow=4, title.position = "top", label.position = "right"))
ST_scale_2_row <- scale_fill_manual(breaks = ST_breaks,values = ST_palette,labels = break_of_sts_renamed, name = "Sequence type (ST)", guide = guide_legend(order=1,nrow=2, title.position = "top", label.position = "right"))

# Sequence type
ST258_color <- subset(ST_palette,names(ST_palette) =="ST258") %>% `names<-`(NULL)
not_ST258_color <- "gray30"
ST258_scale_fill <- scale_fill_manual(breaks = c("ST258","Not ST258"),labels = str_pad(c("ST258","Not ST258"),7,side = 'right'),values = c(ST258_color,not_ST258_color),name = "Sequence type (ST)")

## ST258 with other 
ST258_scale_w_overall_fill <- scale_fill_manual(breaks = c('Overall',"ST258","Not ST258"),labels = str_pad(c("Overall","ST258","Not ST258"),7,side = 'right'),values = c("black",ST258_color,not_ST258_color),name = "Sequence type (ST)")

# Carbapenemase scale
carbapenemase_palette <- readRDS("./lib/carbapenemase_palette.RDS")
carbapenemase_breaks <-  names(carbapenemase_palette)
carbapenemase_colors <-carbapenemase_palette

carbapenemase_scale <- scale_fill_manual(breaks = carbapenemase_breaks,values = carbapenemase_colors,name = "Carbapenemase", guide = guide_legend(ncol=4, title.position = "top", label.position = "right",order=2)) 
carbapenemase_scale_color <- scale_color_manual(breaks = carbapenemase_breaks,values = carbapenemase_colors,name = "Carbapenemase", guide = guide_legend(ncol=4, title.position = "top", label.position = "right",order=2),na.translate = FALSE) 
carbapenemase_scale_4_row <- scale_fill_manual(breaks = carbapenemase_breaks,values = carbapenemase_colors,name = "Carbapenemase", guide = guide_legend(nrow=4, title.position = "top", label.position = "right",order=2)) 
carbapenemase_scale_3_row <- scale_fill_manual(breaks = carbapenemase_breaks,values = carbapenemase_colors,name = "Carbapenemase", guide = guide_legend(nrow=3, title.position = "top", label.position = "right",order=2)) 
carbapenemase_scale_2_row <- scale_fill_manual(breaks = carbapenemase_breaks,values = carbapenemase_colors,name = "Carbapenemase", guide = guide_legend(nrow=2, title.position = "top", label.position = "right",order=2)) 

# Carbapenemase plasmid scale
carbapenemase_plasmid_palette <- readRDS("./lib/carbapenemase_plasmid_palette.RDS")
carbapenemase_plasmid_breaks <-  names(carbapenemase_plasmid_palette)
carbapenemase_plasmid_colors <- carbapenemase_plasmid_palette

plasmid_fill_scale <- scale_fill_manual(breaks = carbapenemase_plasmid_breaks,values = carbapenemase_plasmid_colors, name="Carbapenemase-containing plasmid",guide = guide_legend(title.position = "top", label.position = "right"))
plasmid_color_scale <- scale_color_manual(breaks = carbapenemase_plasmid_breaks,values = carbapenemase_plasmid_colors, name="Carbapenemase-containing plasmid",guide = guide_legend(title.position = "top", label.position = "right"),na.translate = FALSE) 
 
# Carbapenemase replicon type palette
carbapenemase_plasmid_replicon_palette <- readRDS("./lib/plasmid_replicon_palette.RDS")
carbapenemase_plasmid_replicon_breaks <-  names(carbapenemase_plasmid_replicon_palette)
carbapenemase_plasmid_replicon_labels <- gsub(",",", ",carbapenemase_plasmid_replicon_breaks)
carbapenemase_plasmid_replicon_colors <- carbapenemase_plasmid_replicon_palette

plasmid_replicon_type_fill_scale <- scale_fill_manual(breaks = carbapenemase_plasmid_replicon_breaks,values = carbapenemase_plasmid_replicon_colors, labels = carbapenemase_plasmid_replicon_labels, name="Plasmid replicon type(s)",guide = guide_legend(title.position = "top", label.position = "right"))
plasmid_replicon_type_color_scale <- scale_color_manual(breaks = carbapenemase_plasmid_replicon_breaks,values = carbapenemase_plasmid_replicon_colors, labels = carbapenemase_plasmid_replicon_labels, name="Plasmid replicon type(s)",guide = guide_legend(title.position = "top", label.position = "right"),na.translate = FALSE) 

# Carbapenemase content scale
KPC_plasmid_scale <-  scale_fill_manual(breaks = c(T,F),values=c("red","gray"),labels = c("Present","Absent"),name = "Carbapenemase content",guide = guide_legend(nrow=2, title.position = "top", label.position = "right"))
KPC_contig_bool <-  scale_fill_manual(breaks = c(T,F),values=c("red","gray"),labels = c("Present","Absent"),name = "Carbapenemase content",guide = guide_legend(nrow=2, title.position = "top", label.position = "right"))

# Phylogeographic scales
## Decade
decade_palette <- readRDS("./lib/decades_palette.RDS")
decade_breaks <- names(decade_palette)
decade_colors <- decade_palette
decade_scale <- scale_fill_manual(breaks = sort(decade_breaks),values = decade_colors,name = "Decade",guide = guide_legend(nrow=3,order=1),drop=FALSE) 

## Continent
continent_palette <- readRDS("./lib/continent_palette.RDS")
continent_breaks <- names(continent_palette)
continent_colors <- continent_palette
continent_scale <- scale_fill_manual(breaks = sort(continent_breaks),values = continent_colors,name = "Continent",guide = guide_legend(nrow=3,order=2),drop=FALSE)

## California scale
california_scale <- scale_fill_manual(breaks = c(TRUE,FALSE),values = c("red","#CCCCCC"),labels = c("California","Other"),name = "Location",guide = guide_legend(nrow=3,order=3),drop=FALSE)
california_scale_color <- scale_color_manual(breaks = c(TRUE,FALSE),values = c("red","#CCCCCC"),labels = c("California","Other"),name = "Location",guide = guide_legend(nrow=3,order=3),drop=FALSE)

## Study scale
study_scale <- scale_fill_manual(breaks = c("2014-15 study","2021-23 study","Public California genomes","Other public genomes"),values = c("#D55E00","#0072B3","#2B2B2B","white"),name = "Type of assembly",guide = guide_legend(nrow=4,order=4),drop=FALSE)
study_scale_color <- scale_color_manual(breaks = c("2014-15 study","2021-23 study","Public genomes"),values = c("#D55E00","#0072B3","#2B2B2B"),name = "California assemblies",guide = guide_legend(nrow=3,order=4),drop=FALSE)
study_scale_fill <- scale_fill_manual(breaks = c("2014-15 study","2021-23 study","Public genomes"),values = c("#D55E00","#0072B3","#2B2B2B"),name = "California assemblies",guide = guide_legend(nrow=3,order=4),drop=FALSE)
study_comparison_scale_fill <- scale_fill_manual(breaks = c("2014-15","2021-23"),values = c("#D55E00","#0072B3"),name = "Study",guide = guide_legend(nrow=1,order=2),drop=FALSE)

# Mobtyper scales
## Mobility
mobility_scale <- scale_fill_manual(values = c("brown","blue","#C0C0C0"),breaks = c("conjugative",'mobilizable','non-mobilizable'),labels = c("Conjugative","Mobilizable","Non-Mobilizable"),name="Predicted mobility",guide = guide_legend(nrow=2, title.position = "top", label.position = "right"))

## Host range
color_range <- hues::iwanthue(n = 6,hmin = 0,hmax = 75)
host_range_scale <- scale_fill_manual(values = c("darkgray", color_range ),breaks = c("Actinomycetota,Bacillota,Pseudomonadota","Klebsiella","Enterobacteriaceae","Enterobacterales"                            ,"Gammaproteobacteria","Salmonella"),labels =  c("Actinomycetota, Bacillota, Pseudomonadota","Klebsiella","Enterobacteriaceae","Enterobacterales"                            ,"Gammaproteobacteria","Salmonella"),name="Observed host range",guide = guide_legend(nrow=2, title.position = "top", label.position = "right"))

# Resistance colors
resistance_cat_colors <- c("Susceptible" = "#005AB5","Intermediate"="#FFC20A","Resistant" = "#DC3220")
resistance_cat_scale <- scale_fill_manual(breaks = c(names(resistance_cat_colors)),values=resistance_cat_colors,labels = names(resistance_cat_colors), name="Resistance category", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))
resistance_cat_scale_two_title_rows <- scale_fill_manual(breaks = c(names(resistance_cat_colors)),values=resistance_cat_colors,labels = names(resistance_cat_colors), name="Resistance\ncategory", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))

MIC_variables <- c("IMI","MERO","CST","FOS",'blbli','CZA',"IR",'MVB',"PLZ",
                   "CT","DLX","ERV","OMC",'FDC')  %>% rev

resistance_vars <- paste0(MIC_variables,"_dich_num")   