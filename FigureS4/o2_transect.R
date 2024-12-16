### plotting effect of globin expression on consumption and oxygen concentration from surface to center of cluster - figure S4 ###

###history###
# none

### what does this script do?###
#1. Allows user to interactively select files to read .csv output files from simulation, reformats data, and plots

### required directory organization ###
# none

### required filenaming and columnnaming scheme ###
# none

# running R version 4.0.0 (2020-04-24)
# load packages

library(ggplot2)
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
library(tidyr)

# interactively choosing directory and data files
wd <- choose.dir(default = NA, caption = "Select Folder")
s.trans <- file.choose() %>% read.csv()  #SuppFigTransect.csv
setwd(wd)

#rename columns

s.trans1 <- pivot_longer(s.trans, cols = Consumption...no.myoglobin:Oxygen...myoglobin, names_to = "x", values_to = "Value") %>% 
  separate(., col = x, into = c("Amount", "Globin"), sep = "\\...")


red <- "#EC7063"
blue <- "#85C1E9"


g1 <- s.trans1[which(s.trans1$Amount == "Consumption"),]
g1$Globin <- factor(g1$Globin, levels = c("no", "myoglobin"))
g2 <- s.trans1[which(s.trans1$Amount == "Oxygen"),]
g2$Globin <- factor(g2$Globin, levels = c("no", "myoglobin"))

consumption <- ggplot(g1, aes(R, Value, col = Globin)) + geom_line(lwd = 1.5) + geom_smooth(se=FALSE, col = NA) + 
  scale_x_continuous(limits = c(10,20), expand = c(0, 0)) + ylim(c(0, 0.25)) + xlab(expression(paste("Cluster radius from center (", mu,"m",")"))) + 
  ylab(expression(O[2]~consumption~rate)) + 
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.4,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top",
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=0.83) + 
  scale_colour_manual(labels = c("0 mM globin", "0.1 mM globin"), values = c(blue, red))

concentration <- ggplot(g2, aes(R, Value, col = Globin)) + geom_line(lwd = 1.5) + geom_smooth(se=FALSE, col = NA) + 
  scale_x_continuous(limits = c(10,20), expand = c(0, 0)) + scale_y_log10(limit = c(0.00005, 1)) + 
  xlab(expression(paste("Cluster radius from center (", mu,"m",")"))) + ylab(expression(O[2]~concentration~mM~(log[10]))) + 
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.4,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top",
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=0.83) + 
  scale_colour_manual(labels = c("0 mM globin", "0.1 mM globin"), values = c(blue, red))

ggsave(filename = "zoomed_oxygenconsumption_v2.pdf", plot = consumption, width = 15, height = 10, units = "in")

ggsave(filename = "zoomed_oxygenconcentration_v2.pdf", plot = concentration, width = 15, height = 10, units = "in")

