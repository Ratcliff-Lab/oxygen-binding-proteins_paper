### plotting globins relieving growth problems ###

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
low.O2 <- file.choose() %>% read.csv() 
high.O2 <- file.choose() %>% read.csv()
setwd(wd)

# reformating data from wide to long, removing "r" character to numerical element, adding globin concentration column
long.low.O2 <- gather(low.O2, Radius, growth.rate, r10:r70, factor_key=TRUE) %>% mutate(globin = "0.01~mM~globin")
long.low.O2[,2] <- sub("^\\D+", "", long.low.O2[,2])

long.high.O2 <- gather(high.O2, Radius, growth.rate, r10:r70, factor_key=TRUE) %>% mutate(globin = "0.1~mM~globin")
long.high.O2[,2] <- sub("^\\D+", "", long.high.O2[,2])

data <- rbind(long.low.O2, long.high.O2) %>% filter(Radius <= 25)

purple <- "#7678ed"
leaf <- "#81B29A"
beige <- "#F2CC8F"
terra <- "#E07A5F"
denim <- "#457b9d"

gr <- ggplot(data, aes(x = Oxygen, y = growth.rate, col = Radius)) + geom_line(lwd = 1.5) +
  geom_smooth(se=FALSE, col = NA) + facet_wrap2(~globin, labeller = label_parsed, axes = "all", remove_labels = "all") + 
  theme(strip.background = element_blank(), 
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top", legend.background = 
          element_rect(fill="transparent", size=0.5),
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 0), panel.spacing.x = unit(2.5, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio = 0.5) +
  scale_color_manual(name = "Radius", labels = c(expression(paste("10 ",mu,"m")), expression(paste("15 ",mu,"m")),expression(paste("20 ",mu,"m")), expression(paste("25 ",mu,"m"))), 
                     values = c(leaf, beige, terra, purple)) +
  scale_x_continuous(limits = c(0, 0.25), breaks = seq(0.05, 0.25, by = 0.1), expand = c(0, 0)) +
  xlab("Oxygen (mM)") + ylab("Fraction of size-based \ngrowth decline alleviated \nby globin expression")


ggsave(filename = "growthrate_globin_v3.pdf", plot = gr, width = 15, height = 10, units = "in")

### plotting various globin concentrations instead; nov6th ###
gp <- file.choose() %>% read.csv()
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
