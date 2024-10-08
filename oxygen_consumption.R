### plotting Tony's model ###

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
smallwd <- choose.dir(default = NA, caption = "Select Folder")
largewd <- choose.dir(default = NA, caption = "Select Folder")  
setwd(wd)

# importing data  
small_file_list <- list.files(smallwd, pattern = ".csv", full.names = TRUE, recursive = FALSE)
small_data <- lapply(small_file_list, read.csv, stringsAsFactors = FALSE) %>% bind_rows() %>% data.frame()
large_file_list <- list.files(largewd, pattern = ".csv", full.names = TRUE, recursive = FALSE)
large_data <- lapply(large_file_list, read.csv, stringsAsFactors = FALSE) %>% bind_rows() %>% data.frame()
all_data <- rbind(large_data, small_data)
all_data$Oxygen.1 <- factor(all_data$Oxygen.1, levels = c("Low", "High"))
all_data$Globin <- factor(all_data$Globin, levels = c("No", "Yes"))

red <- "#EC7063"
blue <- "#85C1E9"

#plotting data
pp <- ggplot(all_data, aes(X, Normalized.Consumption, col = Oxygen.1, linetype = Globin)) + geom_line(lwd = 1.5) + geom_smooth(aes(group = Oxygen.1),se=FALSE, col = NA)+ 
  facet_grid(~Size) + xlab(expression(paste("Cluster radius from center (", mu,"m",")"))) + ylab(expression(Normalized ~O[2]~ consumption)) + 
  theme(strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "horizontal", legend.justification = c(0.75,1), legend.title = element_blank(), 
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), legend.position = "top",
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio=1) + 
  ylim(c(0,1)) + scale_colour_manual(labels = c(expression(Intermediate ~ O[2]), expression(Supplemental ~ O[2])), values = c(blue, red)) +
  scale_linetype_manual(labels = c("No globin", "Globin"), values = c('dashed','solid'))

#saving plot
ggsave(filename = "oxygen_consumption_v1.pdf", plot = pp, width = 15, height = 10, units = "in")


####### plotting Tony's new data, Apr. 8th, 2023 #########
# high oxygen transect
wd <- choose.dir(default = NA, caption = "Select Folder")
setwd(wd)
O25_G0.1_R_20 <- read.csv("Tonys_modeldata/new_model_8apr2023/O25_G_0.1_R_20.csv")
o.df <- melt(O25_G0.1_R_20 ,  id.vars = 'R', variable.name = 'Groups')

blue.2 <- "#419EDC"
green <- "#A7C270"
red.2 <- "#E95949"
purple <- "#9C80C6"

highO2 <- ggplot(o.df, aes(R, value)) + geom_line(aes(colour = Groups), lwd = 1.5) + xlab(expression(paste("Cluster radius from center (", mu,"m",")"))) + 
  ylab(expression(atop(~O[2]~ consumption ~(mM~L^"-1"~s^"-1")~ and, ~concentration ~(mM~L^"-1")))) +
  theme(strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.75,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top",
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio=1) + 
  scale_colour_manual(labels = c(expression(O[2]~ consumption), expression(O[2]~ concentration), "Unbound myoglobin", expression(O[2]~ bound~myoglobin)), values = c(red.2, blue.2, green, purple)) +
  guides(color = guide_legend(ncol = 2)) 

# 10 micron and 20 micron
small.10 <- file.choose() %>% read.csv() 
  #read.csv("Tonys_modeldata/new_model_8apr2023/new_model_8apr2023/10_micron.csv")
large.20 <- file.choose() %>% read.csv() 
  #read.csv("Tonys_modeldata/new_model_8apr2023/new_model_8apr2023/20_micron.csv")

s.df <- pivot_longer(small.10, 2:5) %>% as.data.frame()
oxygen <- c(rep("High", 44), rep("Low", 44))
globin <- rbind(rep(c(rep("0.1", 22), rep("0", 22)), times = 2)) %>% as.data.frame() %>% t()
s.df <- cbind(s.df[order(s.df$name, decreasing = F),], oxygen, globin, size = rep("Small", nrow(s.df)))

#s <- ggplot(s.df, aes(R, value, col = oxygen, linetype = globin)) + geom_line(lwd = 1.5) + geom_smooth(aes(group = oxygen),se=FALSE, col = NA)


l.df <- pivot_longer(large.20, 2:5) %>% as.data.frame()
oxygen <- c(rep("High", 84), rep("Low", 84))
globin <- rbind(rep(c(rep("0.1", 42), rep("0", 42)), times = 2)) %>% as.data.frame() %>% t()
l.df <- cbind(l.df[order(l.df$name, decreasing = F),], oxygen, globin, size = rep("Large", nrow(l.df)))

#l <- ggplot(l.df, aes(R, value, col = oxygen.l, linetype = globin.l)) + geom_line(lwd = 1.5) + geom_smooth(aes(group = oxygen.l),se=FALSE, col = NA)

f <- rbind (s.df, l.df)
f$normalized <- (f$value)/0.237
f$size <- as.factor(f$size)
f$oxygen <- as.factor(f$oxygen)
f$oxygen <- factor(f$oxygen, levels = c("Low", "High"))
#f$size <- factor(f$size, levels = c("Large", "Small"))
#levels(f$size) <- c("Large", "Small")

red <- "#EC7063"
blue <- "#85C1E9"

s.l.O2 <- ggplot(f, aes(R, normalized, col = oxygen, linetype = globin)) + geom_line(lwd = 1.5) + geom_smooth(aes(group = oxygen),se=FALSE, col = NA)+ 
  facet_grid(~size) + xlab(expression(paste("Cluster radius from center (", mu,"m",")"))) + ylab(expression(atop(Fraction~of~maximum~O[2],~consumption~rate~(mM~L^"-1"~s^"-1")))) + 
  theme(strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.4,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top",
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio=1) + 
  scale_colour_manual(labels = c(expression(0.01~mM ~ O[2]), expression(0.25~mM ~ O[2])), values = c(blue, red)) +
  scale_linetype_manual(labels = c("0 mM myoglobin", "0.1 mM myoglobin"), values = c('dashed','solid')) + ylim(c(0,1)) 

ggsave(filename = "highoxygen_consumption_v1.pdf", plot = highO2, width = 15, height = 10, units = "in")
ggsave(filename = "oxygen_consumption_v4.pdf", plot = s.l.O2, width = 15, height = 10, units = "in")

