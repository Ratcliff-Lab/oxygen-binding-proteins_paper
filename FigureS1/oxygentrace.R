### plotting oxygen saturation curves ###

###history###
# none

### what does this script do?###
#1. Allows user to interactively select .csv files to read and plot

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

# interactively choosing directory and data files
wd <- choose.dir(default = NA, caption = "Select Folder")
o2 <- file.choose() %>% read.csv() # SuppFigOxTraces.csv
setwd(wd)

long.o2 <- o2[,c(1, 12:13)]
colnames(long.o2) <- c("Time", "Low", "Supplemental")
long.o2 <- gather(long.o2, o2, measurement, Low:Supplemental, factor_key=TRUE)

red <- "#EC7063"
blue <- "#85C1E9"

o2.trace <- ggplot(long.o2, aes(x = Time, y = measurement, col = o2)) + geom_line(lwd = 1.5) + xlab("Time (hrs)") +
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.7,1), legend.title=element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        legend.position= "top", legend.background = element_rect(fill="white"),
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 0), panel.spacing.x = unit(1, "lines"),
        axis.ticks=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect("transparent"),
        aspect.ratio = 1) +
  ylab(expression(Percent~saturation~of ~ O[2])) + scale_color_manual(values = c(blue, red), 
    labels = c(expression(No~supplemental ~ O[2]), expression(Supplemental ~ O[2])))

ggsave(filename = "oxygenlevels_v1.pdf", plot = o2.trace, width = 15, height = 10, units = "in")
