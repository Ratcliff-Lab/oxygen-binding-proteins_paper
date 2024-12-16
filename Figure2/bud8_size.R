### plotting size distribution of wildtype and bud8 clusters ###

###history###
# none

### what does this script do?###
#1. Allows user to interactively select files to read .csv files
#2. Uses raw input data to calculate population distribution of particle sizes, using a sliding window
#3. Plots figure 2A

### required directory organization ###
# None, chooses specific files

### required filenaming and columnnaming scheme ###
# None

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
theme_set(theme_pubr())
library(cowplot)
library(rstatix)
library(tidyr)
library(ggprism)
library(dgof)
library(mixtools)

# set working directory 
wd <- choose.dir(default = NA, caption = "Select Folder")
wd1 <- choose.dir(default = NA, caption = "Select Folder For Plots")
setwd(wd1)

# import datafile
bud8.data <- file.choose() %>% read.csv() # bud8_size.csv

ggplot(bud8.data, aes(Radius)) + geom_line(aes(y = Ace2.fraction, col = "red")) + 
  geom_line(aes(y = Ace2.bud8.fraction, col = "blue")) + scale_x_log10(limits = c(1,100))

w = 4
Ace2.smooth <- rep(0, nrow(bud8.data))
Ace2.bud8.smooth <- rep(0, nrow(bud8.data))

for (i in (w+1):(nrow(bud8.data)-w)){
  Ace2.smooth[[i]] <- mean(bud8.data$Ace2.fraction[(i-w):(i+w)])
  Ace2.bud8.smooth[[i]] <- mean(bud8.data$Ace2.bud8.fraction[(i-w):(i+w)])
}

bud8.data <- cbind(bud8.data, Ace2.smooth, Ace2.bud8.smooth)

bs <- ggplot(bud8.data, aes(Radius)) + geom_line(aes(y = Ace2.smooth, col = "#3fa4e8"), lwd = 1.5) + 
  geom_line(aes(y = Ace2.bud8.smooth, col = "#ed1c1c"), lwd = 1.5) + xlab(expression(paste("Radius ","(",mu, "m)"))) + 
  ylab("Frequency") + scale_color_manual(values = c("#3fa4e8", "#ed1c1c"), labels = c(bquote(atop("Normal snowflake","("*italic("ace2")*Delta*")")), bquote(atop("Small snowflake","("*italic("ace2")*Delta~"+"~italic("bud8")*Delta*")")))) +
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.9,1), legend.title=element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        legend.position= "top",
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 0), panel.spacing.x = unit(2, "lines"),
        axis.ticks=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA), aspect.ratio = 1)

ggsave(filename = "small_normalsize_v3.pdf", plot = bs, width = 12, height = 10, units = "in")


