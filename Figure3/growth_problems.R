### plotting globins relieving growth problems - figure 3E ###

###history###
# none

### what does this script do?###
#1. Allows user to interactively select files to read analyzed .csv files from modeling
#2. Plots figure

### required directory organization ###
# none

### required filenaming and columnnaming scheme ###
# none

# running R version 4.0.0 (2020-04-24)
# load packages

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(remotes)
library(shiny)
library(shinyDirectoryInput)
library(magrittr)
library(ggh4x)
library(ggridges)
library(ggpubr)
library(reshape2)
library(cowplot)

### plotting various globin concentrations instead; nov6th ###
gp <- file.choose() %>% read.csv() # nov6_2023.csv
  #"/Tonys_modeldata/growth_problems/nov6_2023.csv"
colnames(gp)[2:6] <- c("0.01", "0.05", "0.10", "0.15", "0.20")
long.gp <- melt(gp, id.vars = "Oxygen", variable.name="Globin", value.name="growth.rate")

purple <- "#7678ed"
leaf <- "#81B29A"
beige <- "#F2CC8F"
terra <- "#E07A5F"
denim <- "#457b9d"

g.r <- ggplot(long.gp, aes(x = Oxygen, y = growth.rate, col = Globin)) + geom_line(lwd = 1.5) +
  geom_smooth(se=FALSE, col = NA) +
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.65,1), legend.title=element_text(size = 34), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1.3,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        legend.position= "top", legend.background = element_rect(fill="white"),
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 0), panel.spacing.x = unit(1, "lines"),
        axis.ticks=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect("transparent"),
        aspect.ratio = 0.3) +
  ylab("Fraction of size-based \ngrowth decline alleviated \nby globin expression") + xlab("Oxygen (mM)") + 
  scale_color_manual(name = "Myoglobin", labels = c("0.01 mM", "0.05 mM", "0.1 mM", "0.15 mM", "0.2 mM"), values = c(leaf, beige, terra, purple, denim)) + 
  scale_x_continuous(expand = c(0, 0)) + guides(color=guide_legend(nrow=2, byrow=TRUE))

ggsave(filename = "growthrate_globin_v6.pdf", plot = g.r, width = 15, height = 10, units = "in")
