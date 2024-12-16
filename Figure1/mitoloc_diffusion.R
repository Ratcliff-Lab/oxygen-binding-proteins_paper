### plotting mitoloc data, testing significance ###

###history###
# none

### what does this script do?###
#1. Allows user to interactively select files to read, analyze for significance, and plot figure
#2. Generates box plot with significance notation

### required directory organization ###
# none

### required filenaming and columnnaming scheme ###
# none

### referenece links for help ###
#https://cran.r-project.org/web/packages/ggprism/vignettes/pvalues.html

# running R version 4.0.0 (2020-04-24)
# load packages

library(ggplot2)
library(dplyr)
library(stringr)
library(remotes)
library(shiny)
library(shinyDirectoryInput)
library(ggpmisc)
library(knitr)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(ggprism)

wd <- choose.dir(default = NA, caption = "Select Folder")
setwd(wd)

#chose the .csv file
ML_low_data_file <- file.choose() %>% read.csv() # Low_Total_notcoded.csv
ML_low_data_file$Group <- factor(ML_low_data_file$Group, levels = c("GOB8", "MyoH", "MyoG"))

ML_supp_data_file <- file.choose() %>% read.csv() # Supplemental_Total_notcoded.csv 
ML_supp_data_file$Group <- factor(ML_supp_data_file$Group, levels = c("GOB8", "MyoH", "MyoG"))

yellow <- "#F1C40F"
red <- "#E74C3C"

#merge df
ML_total <- rbind(ML_low_data_file,ML_supp_data_file)
ML_total$Treatment <- as.factor(ML_total$Treatment)
#levels(ML_total$Treatment) <- c("Intermediate[2]", "Supplemental[2]")
ML_total$Treatment <- factor(ML_total$Treatment, labels = c("Low~O[2]", "Supplemental~O[2]"))

# both: final plot with post-hoc notation
xt <- ML_total %>% rstatix::group_by(Treatment) %>% rstatix::tukey_hsd(Diffusion ~ Group) %>% rstatix::add_xy_position()
  # ML_total %>% group_by(Treatment) %>% anova_test(Diffusion ~ Group)

###
ML_f <- ggplot(ML_total, aes(x = Group, y = Diffusion)) + geom_boxplot(aes(fill = Group), outlier.shape = NA, lwd = 1) + xlab(~O[2]~ "Level") + 
  geom_jitter(color="black",size=1.5,position = position_jitter(width = .2)) +
  ylab(expression("Diffusion depth ("*mu*m*")")) + facet_grid(~Treatment, switch = "both") + 
  add_pvalue(xt[1:2,], label = "p.adj.signif", label.size = 16, tip.length = 0, bracket.size = 1.2, step.increase = 0.05) +
  theme(axis.text.x=element_blank()) + scale_fill_manual(labels = c("GOB8" = "Wild type", "MyoH" = "Myohemerythrin", "MyoG" = "Myoglobin"), 
                                                         values = c("white", yellow, red)) +
  theme(strip.background = element_blank(), legend.key = element_rect(fill = NA, color = NA),
        strip.placement = "outside", strip.text = element_text(size=34),
        legend.key.size = unit(1.5, 'cm'), plot.title = element_text(size = 30),
        axis.title = element_text(size=34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")),
        legend.title = element_blank(), axis.text=element_text(size=34, color = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"), panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 1), axis.text.x = element_blank(), axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y = element_text(vjust = 4), axis.ticks.y=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        panel.spacing.x = unit(1, "lines"),legend.position = "top", legend.justification = c(0.75,1),
        plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=2) + 
  facet_grid(~Treatment, switch = "both", labeller = label_parsed)

ML_f

ggsave(filename = "mitoloc_diffusion_v4.pdf", plot = ML_f, width = 12, height = 10, units = "in")


### extra
# finding max diffusion datapoint in MyoH, supplemental
#####max(subset(ML_total, Group=="GOB8" & Treatment=="Supplemental~O[2]")$Diffusion)

