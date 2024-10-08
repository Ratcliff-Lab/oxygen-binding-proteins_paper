### plotting relative fitness ###

###history###
# from 'relative_fitness_box_v2.R'# script
# Version: v3

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
anc_low_O2_relative_file <- file.choose() %>% read.csv() 
anc_supp_O2_relative_file <- file.choose() %>% read.csv()
anc_df <- rbind(anc_low_O2_relative_file, anc_supp_O2_relative_file)
anc_df$cluster.size <- rep("Normal SF", length(anc_df$Competing.Group))
#anc_df$cluster.size <- rep(expression(paste("Normal Snowflake (ace2",delta,")")), length(anc_df$Competing.Group))

small_low_O2_relative_file <- file.choose() %>% read.csv()
small_supp_O2_relative_file <- file.choose() %>% read.csv()
small_df <- rbind(small_low_O2_relative_file, small_supp_O2_relative_file)
small_df$cluster.size <- rep("Small SF", length(small_df$Competing.Group))
#small_df$cluster.size <- rep(expression("Small Snowflake (ace2\delta + bud8\delta)"), length(small_df$Competing.Group))

unicell_lowO2_relative_file <- file.choose() %>% read.csv()
unicell_suppO2_relative_file <- file.choose() %>% read.csv()
unicell_df <- rbind(unicell_lowO2_relative_file, unicell_suppO2_relative_file)
unicell_df$Relative.Fitness <- unicell_df$Relative.Fitness - 1
unicell_df$cluster.size <- rep("Unicellular yeast", length(unicell_df$Competing.Group))

fitness_df <- rbind(anc_df, small_df, unicell_df)
fitness_df$cluster.size <- as.factor(fitness_df$cluster.size)
#fitness_df$cluster.size <- factor(fitness_df$cluster.size, labels = c("Normal Snowflake (ace2~delta)", "Small Snowflake (ace2~delta + bud8~delta)"))
## renaming cluster size factor
fitness_df$cluster.size <- factor(fitness_df$cluster.size, labels = c(bquote(atop("Normal snowflake","("*italic("ace2")*Delta*")")), bquote(atop("Small snowflake","("*italic("ace2")*Delta~"+"~italic("bud8")*Delta*")")), bquote(atop("Unicellular yeast"))))
fitness_df$Relative.Fitness <- fitness_df$Relative.Fitness + 1

fitness_df$Treatment <- as.factor(fitness_df$Treatment)
fitness_df$Competing.Group <- as.factor(fitness_df$Competing.Group)

### dot and SD version
## first find the mean and SD for each group
fitness_df_stat <- fitness_df %>%                   
  group_by(Competing.Group, Treatment, cluster.size) %>%
  summarise_at(vars(Relative.Fitness),list(mean = mean, sd = sd)) %>% as.data.frame()

### paired t-test
## regroup data, merge strain and size
y <- fitness_df %>% group_by(cluster.size, Competing.Group, Treatment) %>% 
  unite('Merged', cluster.size,Competing.Group) %>% as.data.frame()
## group data and analyze
p.paired <- y %>% group_by(Merged) %>% t_test(data =., Relative.Fitness ~ Treatment) %>% add_significance("p") %>% add_y_position()
# ungroup column
p.paired <- separate(p.paired, col = Merged, into = c("Size", "Genotype"), sep = "\\_")
p.paired <- rename(p.paired, Competing.Group = Genotype, cluster.size = Size) 
p.paired$cluster.size <- as.factor(p.paired$cluster.size) #
p.paired$Competing.Group <- as.factor(p.paired$Competing.Group) #
p.paired$group1 <- c("GFP2 - myoH", "GFP2+HEM - myoG+HEM", "GFP2 - myoH", "GFP2+HEM - myoG+HEM", "GFP2 - myoH", "GFP2 - Y55", "GFP2+HEM - myoG+HEM")
p.paired$group2 <- c("GFP2 - myoH", "GFP2+HEM - myoG+HEM", "GFP2 - myoH", "GFP2+HEM - myoG+HEM", "GFP2 - myoH", "GFP2 - Y55", "GFP2+HEM - myoG+HEM")
# adjusting position
p.paired$y.position <- p.paired$y.position + 0.001

### One-sample t-test
x <- fitness_df %>%                   
  group_by(Competing.Group, Treatment, cluster.size) %>% t_test(Relative.Fitness ~ 1, mu = 1) %>% add_significance("p")
x$Competing.Group <- factor(x$Competing.Group, levels = c("GFP2 - myoH", "GFP2+HEM - myoG+HEM", "GFP2 - Y55"))
x1 <- x[order(x$Competing.Group),]


# table
# df <- data.frame("O2 carrier" = rep(c("Myohemerytherin", "Myoglobin"), each = 4), 
#                  "O2 level" = rep(c("Intermediate", "Intermediate", "Supplemental","Supplemental"), times = 2),
#                  "Cluster size" = rep(c("Normal", "Small"), times = 4), "p-value" = x$p, "statistic" = x$statistic, "p.signif" = x$p.signif)
# dt <- ggtexttable(df, rows = NULL, theme = ttheme("light")) 

# colors and annotations for plotting
red <- "#EC7063"
blue <- "#85C1E9"
# dataframe for adding significance annotations
rel.fit.annot <- data.frame(x1 = c(0.8, 1.8, 0.6, 1.6, 2.6), x2 = c(1.2, 2.2, 1.4, 2.4, 3.4), 
                            y1 = c(1.03725, 1.03325, 1.009, 1.0125, 1.0065), y2 = c(1.03725,1.03325,1.009,1.0125,1.0065),
                            cluster.size = c('atop("Normal snowflake", "(" * italic("ace2") * Delta * ")")', 
                            'atop("Normal snowflake", "(" * italic("ace2") * Delta * ")")', 
                            'atop("Unicellular yeast")', 'atop("Unicellular yeast")', 'atop("Unicellular yeast")'))



## plot mean and SD range
rel.fit <- 
  ggplot(fitness_df_stat, aes(x=Competing.Group, y=mean, color=Treatment)) + 
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),position = position_dodge(width=0.5), size = 2) +
  #creating a panel of plots
  facet_grid(~cluster.size, scales = "free") + scale_x_discrete(labels = c("Myohemerythrin", "Myoglobin", "Y55")) + 
  stat_pvalue_manual(p.paired[c(1:2,5:7),], label = "p.signif", label.size = 16, tip.length = 0) + 
  # cleaner formatting
  theme(strip.background = element_blank(), legend.key = element_rect(fill = NA),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), plot.title = element_text(size = 30),
        axis.title = element_text(size =34), legend.text = element_text(size=30), axis.title.y = element_text(vjust = 4),
        legend.title = element_text(size=34), axis.text=element_text(size=34, color = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"), panel.background = element_rect("white"),
        axis.line = element_line(color = "black", size = 1), axis.title.x=element_blank(), 
        axis.text.x = element_text(size = 34,angle = 30, vjust = 0.75, hjust=0.75), 
        legend.position= c(0.85,0.8), legend.background = element_rect(fill="white", size=0.5, linetype="solid", color = "black"),
        axis.ticks=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        panel.spacing.x = unit(1, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1.5) + 
  # horizontal line at 1
  geom_hline(yintercept=1, linetype = "dashed") + ylab("Relative fitness during growth") + 
  scale_colour_manual(name = ~O[2]~ "level", labels = c("Low", "Supplemental"), values = c(blue, red)) + 
  facet_grid(~cluster.size, labeller = label_parsed, scales = "free") + ylim(c(0.99,1.04)) + 
  #annotating 
  geom_segment(data = rel.fit.annot, aes(x = x1, xend = x1, y = y1, yend = y2), colour = "black") +
  geom_segment(data = rel.fit.annot, aes(x = x2, xend = x2, y = y1, yend = y2), colour = "black") + 
  geom_segment(data = rel.fit.annot, aes(x = x1, xend = x2, y = y2, yend = y2),colour = "black")


rel.fit

# saving output as .pdf
ggsave(filename = "relative_fit_anc_v6.pdf", plot = rel.fit, width = 15, height = 12, units = "in")
ggsave(filename = "one_sampletest_fitness_v1.pdf", plot = dt, width = 5, height = 4, units = "in")


