### plotting mitoloc data, testing significance ###

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
ML_low_data_file <- file.choose() %>% read.csv()
ML_low_data_file$Group <- factor(ML_low_data_file$Group, levels = c("GOB8", "MyoH", "MyoG"))

ML_supp_data_file <- file.choose() %>% read.csv()
ML_supp_data_file$Group <- factor(ML_supp_data_file$Group, levels = c("GOB8", "MyoH", "MyoG"))

yellow <- "#F1C40F"
red <- "#E74C3C"

#merge df
ML_total <- rbind(ML_low_data_file,ML_supp_data_file)
ML_total$Treatment <- as.factor(ML_total$Treatment)
#levels(ML_total$Treatment) <- c("Intermediate[2]", "Supplemental[2]")
ML_total$Treatment <- factor(ML_total$Treatment, labels = c("Low~O[2]", "Supplemental~O[2]"))

#y <- c(`Low` = "Intermediate", `Supplemental` = "Supplemental O")
#levels(ML_total$Treatment) <- c('"Intermediate "*O[2]', '"Supplemental "*O[2]')

### testing data
# one.way <- aov(Diffusion ~ Group, data = ML_low_data_file)
# #summary(one.way)
# TR <- TukeyHSD(one.way, conf.level = 0.95)
#plot(TukeyHSD(one.way))

#global anova
# oneway <- compare_means(Diffusion ~ Group,  data = ML_low_data_file, method = "anova")
# ggplot(ML_low_data_file, aes(x = Group, y = Diffusion, fill = Group)) + geom_boxplot() + 
#   stat_compare_means(method = "anova")

# using rstatix package
# x <- tukey_hsd(ML_low_data_file, Diffusion ~ Group) %>% add_xy_position()

#comparison to reference
# ML <- ggplot(ML_low_data_file, aes(x = Group, y = Diffusion)) + geom_boxplot(aes(fill = Group)) + 
#   add_pvalue(x, label = "p.adj.signif", y.position = c(37,39,42))
# 
# ML
# 
# xs <- x <- tukey_hsd(ML_supp_data_file, Diffusion ~ Group) %>% add_xy_position()
# 
# ggplot(ML_total, aes(x = Treatment, y = Diffusion)) + geom_boxplot(aes(fill = Group)) + 
#   add_pvalue(x, label = "p.adj.signif", y.position = c(37,39,42)) + add_pvalue(xs, label = "p.adj.signif", y.position = c(37,39,42))

# both: final plot with post-hoc notation
xt <- ML_total %>% rstatix::group_by(Treatment) %>% rstatix::tukey_hsd(Diffusion ~ Group) %>% rstatix::add_xy_position()
  # ML_total %>% group_by(Treatment) %>% anova_test(Diffusion ~ Group)

# remove outliers
# find.limits <- function(x) {
#   upper.limit <- quantile(x)[3] + 1.5*IQR(x)
#   lower.limit <- quantile(x)[3] - 1.5*IQR(x)
#   ranges <- c(lower.limit, upper.limit)
#   return(ranges)
# }
# 
# MM <- ML_total %>% group_by(Group, Treatment) %>% summarise(MM = find.limits(Diffusion))
# 
# UL <- max(MM$MM)
# LL <- min(MM$MM)

### removing outlier points 
# ML_f <- ML_total %>%
#   group_by(Group, Treatment) %>%
#   mutate(Diffusion.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
#     between(Diffusion, 
#             quantile(Diffusion)[2] - 1.5*IQR(Diffusion),
#             quantile(Diffusion)[4] + 1.5*IQR(Diffusion)))) %>% 
#   ggplot(aes(Group, Diffusion)) + 
#   geom_boxplot(aes(fill = Group), outlier.shape = NA) + 
#   geom_jitter(aes(alpha=Diffusion.show), show.legend=FALSE, width = 0.25) +
#   scale_alpha_continuous(range = c(0, 1)) + facet_grid(~Treatment, switch = "both") + 
#   add_pvalue(xt, label = "p.adj.signif", label.size = 8, tip.length = 0, y.position=c(44,47,50)) +
#   ylab(expression("Diffusion Depth,"~mu*m)) + 
#   scale_fill_manual(name = "Strain", labels = c("GOB8" = "GOB8", "MyoH" = "Myohemerythrin", "MyoG" = "Myoglobin"), values = c("white", yellow, red)) + 
#   theme(strip.background = element_blank(), legend.key = element_rect(fill = NA),
#         strip.placement = "outside", strip.text = element_text(size=22, vjust = 0),
#         legend.key.size = unit(1, 'cm'), plot.title = element_text(size = 24),
#         axis.title = element_text(size =22), legend.text = element_text(size=20), axis.title.y = element_text(vjust = 4),
#         legend.title = element_text(size=22), axis.text=element_text(size=22, color = "black"),
#         plot.margin = margin(20, 150, 20, 150), panel.background = element_rect("white"),
#         axis.line = element_line(color = "black"), axis.text.x = element_blank(), axis.title.x=element_blank(),
#         axis.ticks.x=element_blank()) + facet_grid(~Treatment, switch = "both", labeller = label_parsed) 

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

### links
#https://stackoverflow.com/questions/61335789/how-to-exclude-outliers-when-using-geom-boxplot-geom-jitter-in-r
#https://cran.r-project.org/web/packages/ggprism/vignettes/pvalues.html

