### plotting relative fitness of induced multicell competitions - figure 4 ###

###history###
# from 'relative_fitness_box_v3.R'# script
# Version: v1 - renamed

### what does this script do?###
#1. Allows user to interactively select files to read analyzed .csv outputs from segmentation_fitness_analysis_v1.R
#1a. Interactively choose each file for each variable, which has different experimental conditions
#2. Renames metadata (experimental conditions)
#3. Calculates mean of 5 replicate samples, one-sample t-test, adds notations to a dot plot (mean +- sd)

### required directory organization ###
# Only .csv files to be plotted in folder

### required filenaming and columnnaming scheme ###
# None, out put from previous script should be correct

### referenece links for help ###
# https://trinkerrstuff.wordpress.com/2018/03/15/2246/
# https://www.r-bloggers.com/2018/11/adding-different-annotation-to-each-facet-in-ggplot/

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
library(rstatix)
library(tidyr)
library(ggprism)
library(cowplot)

wd <- choose.dir(default = NA, caption = "Select Folder") # choose where you want plot to end up
setwd(wd)

# interactively choose file and read csv file
ind_0_low_file <- file.choose() %>% read.csv() # Analysis0low_20241009.csv
ind_0_supp_file <- file.choose() %>% read.csv() # Analysis0high_20241009.csv
ind_10_low_file <- file.choose() %>% read.csv() # Analysis10low_20241009.csv
ind_10_supp_file <- file.choose() %>% read.csv() # Analysis10high_20241009.csv
ind_20_low_file <- file.choose() %>% read.csv() # Analysis20low_20241009.csv
ind_20_supp_file <- file.choose() %>% read.csv() # Analysis20high_20241009.csv
ind_50_low_file <- file.choose() %>% read.csv() # Analysis50low_20241009.csv
ind_50_supp_file <- file.choose() %>% read.csv() # Analysis50high_20241009.csv
ind_200_low_file <- file.choose() %>% read.csv() # Analysis200low_20241009.csv
ind_200_supp_file <- file.choose() %>% read.csv() # Analysis200high_20241009.csv

ind_mc_data <- rbind(ind_0_low_file, ind_0_supp_file, ind_10_low_file, ind_10_supp_file, ind_20_low_file, ind_20_supp_file,
                     ind_50_low_file, ind_50_supp_file, ind_200_low_file, ind_200_supp_file)

ind_mc_data$Rel.Fitness <- ind_mc_data$Rel.Fitness + 1

ind_mc_data$Treatment <- as.factor(ind_mc_data$Treatment)

### dot and SD version
## first find the mean and SD for each group
fitness_df_stat <- ind_mc_data %>%                   
  group_by(Induction, Treatment) %>%
  summarise_at(vars(Rel.Fitness),list(mean = mean, sd = sd)) %>% as.data.frame()
fitness_df_stat$Induction <- factor(fitness_df_stat$Induction, levels = c("0", "10", "20", "50", "200"))

### paired t-test
## regroup data, merge strain and size
y <- ind_mc_data %>% group_by(Induction, Treatment) %>% as.data.frame() # before %>% unite('Merged', cluster.size,Competing.Group) %>% as.data.frame()
## group data and analyze
p.paired <- y %>% group_by(Induction) %>% t_test(data =., Rel.Fitness ~ Treatment) %>% add_significance("p") %>% add_y_position()
# ungroup column
# p.paired <- separate(p.paired, col = Merged, into = c("Size", "Genotype"), sep = "\\_")
# p.paired <- rename(p.paired, Competing.Group = Genotype, cluster.size = Size) 
p.paired$Induction <- factor(p.paired$Induction, levels = c("0", "10", "20", "50", "200")) #
# p.paired$Competing.Group <- as.factor(p.paired$Competing.Group) #
p.paired$group1 <- c("0", "10", "20", "50", "200")
p.paired$group2 <- c("0", "10", "20", "50", "200")
# adjusting position
p.paired$y.position <- p.paired$y.position + 0.001

### One-sample t-test
x <- ind_mc_data %>%                   
  group_by(Induction, Treatment) %>% t_test(Rel.Fitness ~ 1, mu = 1) %>% add_significance("p")


# table
# df <- data.frame("O2 carrier" = rep(c("Myohemerytherin", "Myoglobin"), each = 4), 
#                  "O2 level" = rep(c("Intermediate", "Intermediate", "Supplemental","Supplemental"), times = 2),
#                  "Cluster size" = rep(c("Normal", "Small"), times = 4), "p-value" = x$p, "statistic" = x$statistic, "p.signif" = x$p.signif)
# dt <- ggtexttable(df, rows = NULL, theme = ttheme("light")) 

# colors and annotations for plotting
red <- "#EC7063"
blue <- "#85C1E9"
# dataframe for adding significance annotations
rel.fit.annot <- data.frame(x1 = c(0.8, 0.7, 0.8), x2 = c(1.2, 1.3, 1.2), 
                            y1 = c(1.008, 1.015, 1.015), y2 = c(1.008,1.015,1.015),
                            Induction = c("10", "20", "50")) 
rel.fit.annot$Induction <- factor(rel.fit.annot$Induction, levels = c("10", "20", "50"))




## plot mean and SD range
rel.fit <- 
  ggplot(fitness_df_stat, aes(x=Induction, y=mean, color=Treatment)) + 
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),position = position_dodge(width=0.5), size = 2) +
  #creating a panel of plots
  facet_grid(~Induction, scales = "free") + scale_x_discrete(labels = c("0", "10", "20", "50", "200")) + 
  stat_pvalue_manual(p.paired[c(2:4),], label = "p.signif", label.size = 16, tip.length = 0) + 
  ylab("Relative fitness during growth") + xlab(expression(paste(mu, "M", " of Estradiol"))) + 
  # cleaner formatting
  theme(strip.background = element_blank(), legend.key = element_rect(fill = NA),
        strip.placement = "bottom", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), plot.title = element_text(size = 30),
        axis.title = element_text(size =34), legend.text = element_text(size=30), axis.title.y = element_text(vjust = 4),
        legend.title = element_text(size=34), axis.text=element_text(size=34, color = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"), panel.background = element_rect("white"),
        axis.line = element_line(color = "black", size = 1), 
        axis.text.x = element_blank(), 
        legend.position= "top", legend.background = element_rect(fill="white", size=0.5, linetype="solid", color = "black"),
        axis.ticks=element_line(size = 1), axis.ticks.length=unit(.25, "cm"), axis.ticks.x = element_blank(),
        panel.spacing.x = unit(1, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=2) + 
  # horizontal line at 1
  geom_hline(yintercept=1, linetype = "dashed") + 
  scale_colour_manual(name = ~O[2]~ "level", labels = c("Low", "Supplemental"), values = c(blue, red)) + 
  facet_grid(~Induction, labeller = label_parsed, scales = "free", switch="both") + ylim(c(0.985,1.02)) + 
  #annotating 
  geom_segment(data = rel.fit.annot, aes(x = x1, xend = x1, y = y1, yend = y2), colour = "black") +
  geom_segment(data = rel.fit.annot, aes(x = x2, xend = x2, y = y1, yend = y2), colour = "black") + 
  geom_segment(data = rel.fit.annot, aes(x = x1, xend = x2, y = y2, yend = y2),colour = "black")


rel.fit 

onesample <- ggtexttable(x, rows = NULL, theme = ttheme("light"))
pairedsample <- ggtexttable(p.paired[c(1, 7:10)], rows = NULL, theme = ttheme("light"))

# saving output as .pdf
ggsave(filename = "induced_MC_fitness_v2.pdf", plot = rel.fit, width = 15, height = 12, units = "in")
ggsave(filename = "one_sampletest_induced_MC_fitness_v1.pdf", plot = onesample, width = 10, height = 4, units = "in")
ggsave(filename = "pairedsampletest_induced_MC_fitness_v1.pdf", plot = pairedsample, width = 10, height = 4, units = "in")
