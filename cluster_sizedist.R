### plotting cluster size data ###

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

# list of full path files
file_list <- list.files(wd, pattern = ".csv", full.names = TRUE, recursive = FALSE)

# reading each file and mergining into one dataframe
all_data <- lapply(file_list, read.csv, stringsAsFactors = FALSE) %>% bind_rows() %>% data.frame()

# renaming strains
all_data$strain[which(all_data$strain == "GFPH")] <- "GFP2 + HEM3"
all_data$strain[which(all_data$strain == "myoH")] <- "Myohemerythrin"
all_data$strain[which(all_data$strain == "myoG")] <- "Myoglobin + HEM3"

# organizing data
all_data$strain <- factor(all_data$strain, levels = c("GFP2", "GFP2 + HEM3", "Myohemerythrin", "Myoglobin + HEM3"))
#all_data$time <- factor(all_data$time, levels = c("t0", "t60", "t120", "t180"))


# adding column 
type <- function(x) {
  ifelse(x == "GFP2" | x == "GFP2 + HEM3", "Control", "O2-Carriers")}

all_data$category <- type(all_data$strain) 
#all_data <- do.call("cbind", category)

all_data$samplenum <- as.character(all_data$samplenum)

# adding back in volume_data column
all_data$volume_data <- sapply(all_data$diameter_data, function(diameter_data) ((diameter_data/2)**3)*pi*(4/3))
list.alldata <- all_data %>% group_split(samplenum, strain, time)

#############################################################################################################################################################################################

#### starting to bin data #####
# log10 of volume data
all_data$volume_data <- log10(all_data$volume_data) ###### important!
# min and max of bin 
all_vol_min <- min(all_data$volume_data)
all_vol_max <- max(all_data$volume_data)
bins <- 100 # of bins created
binwidth <- (all_vol_max - all_vol_min)/bins # binwidth based on number of bins

# sorts datapoints into bins
all_data_bin <- all_data %>% mutate(all_vol_bin = cut(volume_data, breaks = seq(min(all_data$volume_data),max(all_data$volume_data), by = binwidth)))
# splitting dataframe based on strain, replicate population, and timepoint
list_all_data_bin <- all_data_bin %>% group_split(samplenum, strain, time)
# function: for each item in list, sort by bin, find the average volume of all points and frequency within bin
all_data_weight <- lapply(list_all_data_bin, function(x){
  #x1 <- x %>% group_by(all_vol_bin) %>% summarise(n = n(), avg = mean(volume_data), weight = sum(volume_density)) %>% as.data.frame 
  #x1 <- x %>% group_by(all_vol_bin) %>% summarise(n = n(), avg = mean(volume_data)) %>% mutate(weight = (n*avg)/sum(x$volume_density)) %>% as.data.frame 
  ### for each bin, find count, average of all volume in bin, weight is the average of all points * count normalized by total population volume
  x1 <- x %>% group_by(all_vol_bin) %>% summarise(n = n(), avg = mean(volume_data), weight = (avg*n)/sum(x$volume_data)) %>% as.data.frame
  meta <- x[1,c(3:6)] %>% as.data.frame
  x1 <- cbind(x1, meta)
})
# convert list output into a dataframe
weighted_data.plot <- do.call("rbind", all_data_weight)
weighted_data.plot$samplenum[which(weighted_data.plot$time == "t0")] <- "Unevolved" #
weighted_data.plot$samplenum <- as.character(weighted_data.plot$samplenum)
weighted_data.plot$time <- factor(weighted_data.plot$time, levels = c("t0", "t60", "t120", "t180")) #

# finding the population mean and variance
z <- lapply(list.alldata, function(x){
  avg <- x %>% summarise(pop.avg = mean(volume_data), pop.var = var(volume_data), pop.coe = (sd(volume_data)/mean(volume_data))) %>% as.data.frame 
  meta <- x[1,c(3:6)] %>% as.data.frame
  avg <- cbind(avg, meta)})
z1 <- do.call("rbind", z)
z1$samplenum[which(z1$time == "t0")] <- "Unevolved"
z1$samplenum <- as.character(z1$samplenum)
#z1$time <- factor(z1$time, levels = c("t0", "t60", "t120", "t180"))
z1$category <- factor(z1$category, levels = c("Control", "O2-Carriers"))
z1$log.pop.avg <- log10(z1$pop.avg)
z1$log.pop.var <- log10(z1$pop.var)


breaks <- 10**c(2:6) # used to be 1:5
# scientific_10 <- function(x) {
#   parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
# }

### plotting every distribution
d3 <- ggplot(weighted_data.plot, aes(x = avg, y = weight, col = samplenum)) + geom_line(position = "identity") + 
  facet_nested(time~category + strain) + theme_classic() + scale_fill_discrete(name="Replicate Populations") + 
  scale_color_discrete(name="Replicate Populations") + ggtitle("180-Day Evolved Populations") + 
  ylab("% of Total Population") + scale_x_continuous(labels  = function(x) format(breaks, scientific = F)) + 
  xlab(expression(paste("Cluster Volume, ", mu,"m"^{3}))) + guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(plot.title = element_text(size = 22),
        axis.title = element_text(size =18), legend.text = element_text(size=12), legend.title = element_text(size=15), 
        axis.text=element_text(size=8), strip.text = element_text(size = 16),
        panel.background = element_rect("white"), axis.line = element_line(color = "black")) 

d3

#### plotting quadrant of ancestor and t180 strains ###
# creating a df of only ancestor and t180

# anc.t180.df <- weighted_data.plot %>% filter(time == "t0" | time == "t180")
anc.t180.df <- weighted_data.plot %>% filter(time == "0" | time == "180")
anc.t180.df$tentoavg <- 10^(anc.t180.df$avg)

# line colors
terra <- "#E07A5F"
denim <- "#457b9d"
leaf <- "#81B29A"
beige <- "#F2CC8F"
purple <- "#7678ed"

## plotting quadrant of t180 cluster distributions

cluster.df <- ggplot(anc.t180.df, aes(x = avg, y = weight, col = samplenum)) + geom_line(position = "identity", lwd = 1) + 
  facet_wrap(vars(strain), nrow = 2) + scale_x_continuous(labels = c(expression("10"^2),expression("10"^3),expression("10"^4),expression("10"^5),expression("10"^6))) + 
  #scale_x_continuous(label = scientific_10(breaks)) +
  #scale_x_continuous(labels  = function(x) format(breaks, scientific = T)) + 
  ylab("Frequency") + xlab(expression(paste("Cluster Volume (", mu,"m"^{3},")"))) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(strip.background = element_blank(), legend.key = element_blank(),
      strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
      legend.key.size = unit(1, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34), 
      axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
      axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
      panel.background = element_rect("white"), legend.position = "none", 
      axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 30), panel.spacing.x = unit(2, "lines"), 
      panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio = 0.75) + 
  scale_color_manual(name = "Replicate Populations", values = c("black", terra, denim, leaf, beige, purple)) 

cluster.df 

ggsave(filename = "evolved_clusterdist_v3.pdf", plot = cluster.df, width = 15, height = 10, units = "in")


### plotting population mean and variance 
labels <- c("0", "60", "120", "180")

pop.avg <- ggplot(z1, aes(time, pop.avg, col = samplenum)) + geom_point(size = 6) + facet_wrap(vars(strain), nrow = 2) + 
  scale_x_discrete(labels = labels) + ylab(expression(paste("Average cluster volume (", mu,"m"^{3},")"))) +
  xlab("Days of evolution") + 
  theme(strip.background = element_blank(), legend.key = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34), 
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), legend.position = "none", 
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_color_manual(name = "Replicate Populations", values = c(terra, denim, leaf, beige, purple, "black")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

pop.var <- ggplot(z1, aes(time, pop.var, col = samplenum)) + geom_point(size = 6) + facet_wrap(vars(strain), nrow = 2) + 
  scale_x_discrete(labels = labels) + ylab(expression(paste("Variance in cluster volume (", mu,"m"^{3},")"))) +
  xlab("Days of evolution") + 
  theme(strip.background = element_blank(), legend.key = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34), 
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), legend.position = "none", 
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_color_manual(name = "Replicate Populations", values = c(terra, denim, leaf, beige, purple, "black")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

pop.coe <- ggplot(z1, aes(time, pop.coe, col = samplenum)) + geom_point(size = 6) + facet_wrap(vars(strain), nrow = 2) + 
  scale_x_discrete(labels = labels) + ylab(expression(paste("Coefficient of variation in cluster volume (", mu,"m"^{3},")"))) +
  xlab("Days of evolution") + 
  theme(strip.background = element_blank(), legend.key = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34), 
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), legend.position = "none", 
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA),
        axis.ticks.y=element_line(size = 1), axis.ticks.x=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        prism.ticks.length.y = unit(0.15, "cm")) + 
  scale_color_manual(name = "Replicate Populations", values = c(terra, denim, leaf, beige, purple, "black")) + 
  scale_y_continuous(guide = guide_prism_minor(), limits = c(0, 3), expand = c(0, 0))

### two-way anova on CoE
CV <- z1 %>% anova_test(pop.coe ~ strain*time) %>% as.data.frame()
CV.tab <- ggtexttable(CV, rows = NULL, theme = ttheme("light"))
ggsave(filename = "CVtab.pdf", plot = CV.tab, width = 5, height = 4, units = "in")


ggsave(filename = "evolved_clusterpopmean_v1.pdf", plot = pop.avg, width = 15, height = 10, units = "in")
ggsave(filename = "evolved_clusterpopvariance_v1.pdf", plot = pop.var, width = 15, height = 10, units = "in")
ggsave(filename = "evolved_clusterpopcoe_v1.pdf", plot = pop.coe, width = 15, height = 10, units = "in")

###### fold-change ######
median_all_data <- lapply(all_data_weight, function(x){
  med <- x %>% filter(cumsum(weight) <= 0.5 | cumsum(weight) < 0.52) %>% summarise(y_med = last(weight), med_bin = last(all_vol_bin), med_Sum = sum(weight))
  nine_five <- x %>% filter(cumsum(weight) <= 0.95 | cumsum(weight) < 0.96) %>% summarise(y_95 = last(weight), ninefive_bin = last(all_vol_bin), ninefive_Sum = sum(weight))
  five <- x %>% filter(cumsum(weight) <= 0.05 | cumsum(weight) < 0.07) %>% summarise(y_05 = last(weight), five_bin = last(all_vol_bin), five_Sum = sum(weight))
  
  meta <- x[1,c(5:8)]
  
  med_bin <- as.character(med[,2])
  nine_five_bin <- as.character(nine_five[,2])
  five_bin <- as.character(five[,2])
  
  med_index <- which(x == med_bin)
  nine_five_index <- which(x == nine_five_bin)
  five_index <- which(x == five_bin)
  
  med_cluster <- ((6/pi)*(10**(x[med_index,3])))**(1/3)
  nine_five_cluster <- ((6/pi)*(10**(x[nine_five_index,3])))**(1/3)
  five_cluster <- ((6/pi)*(10**(x[five_index,3])))**(1/3)
  
  perc <- cbind(meta, med, med_index, med_cluster, nine_five, nine_five_index, nine_five_cluster, five, five_index, five_cluster)
})

median_all_data.plot <- do.call("rbind", median_all_data)

median_all_data.plot$time <- factor(median_all_data.plot$time, levels = c("t0", "t60", "t120", "t180"))
median_all_data.plot$strain <- factor(median_all_data.plot$strain, levels = c("GFP2", "GFPH", "myoH", "myoG")) ###

# med.clust <- ggplot(median_all_data.plot[,c(1:9)], aes(x = time, y = med_cluster, shape = samplenum, group = strain)) + geom_point(aes(col = strain), position=position_dodge(width=0.5), size = 5) + 
#   scale_y_continuous(breaks = seq(0,60,by=15), limits=c(0,60)) + theme_classic() + theme(panel.grid.major.y = element_line(color = "gray",size = 0.5,linetype = 2), 
#                                                                                          axis.title = element_text(size =15)) + xlab("Time") + ylab(expression(paste("Cluster Diameter, ", mu,"m"))) + scale_color_discrete(name="Strains") + 
#   scale_shape_discrete(name = "Replicate Populations") + ggtitle("50th-Percentile Cluster Size")

#ggplot(median_all_data.plot, aes(x = samplenum, y = med_cluster, col = samplenum, group = time)) + geom_point(size = 5) + geom_errorbar(aes(ymin=five_cluster, ymax=nine_five_cluster), width = 0.5) + facet_nested(~strain+time, switch = "x")

### 50th percentile ### 
fold <- lapply(1:nrow(median_all_data.plot), FUN = function(i) {
  s <- as.character(median_all_data.plot$strain[i])
  unevolved <- median_all_data.plot$med_cluster[which(median_all_data.plot$strain == s & median_all_data.plot$samplenum == 0)]
  fold.change <- c(s, median_all_data.plot$samplenum[i], as.character(median_all_data.plot$time[i]), median_all_data.plot$category[i], (median_all_data.plot$med_cluster[i]/unevolved)-1)
})

fifty.fold.change.data <- data.frame(do.call(rbind, fold), stringsAsFactors = F)
colnames(fifty.fold.change.data) <- append(colnames(median_all_data.plot[1:4]), "fold.data")
fifty.fold.change.data$fold.data <- as.numeric(fifty.fold.change.data$fold.data)

fifty.fold.change.data$time <- factor(fifty.fold.change.data$time, levels = c("t0", "t60", "t120", "t180"))
fifty.fold.change.data$strain <- factor(fifty.fold.change.data$strain, levels = c("GFP2", "GFPH", "myoH", "myoG")) ###

d4 <- ggplot(fifty.fold.change.data, aes(time, fold.data)) + geom_boxplot(aes(col = strain), outlier.shape = NA) + 
  geom_jitter(width = 0.1) + geom_hline(yintercept = 0, linetype = "dotted") + facet_wrap(~strain, ncol = 2) + 
  theme(plot.title = element_text(size = 22), legend.position='right', strip.background = element_blank(),
        strip.text.x = element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=18), axis.title = element_text(size =18)) + 
  ylab("Fold change in cluster diameter") + xlab("Time") + scale_color_discrete(name="Strains") + ggtitle("50th-Percentile Cluster Diameter of Evolved Populations")

d4

### 95th percentile ###
fold <- lapply(1:nrow(median_all_data.plot), FUN = function(i) {
  s <- as.character(median_all_data.plot$strain[i])
  unevolved <- median_all_data.plot$nine_five_cluster[which(median_all_data.plot$strain == s & median_all_data.plot$samplenum == 0)]
  fold.change <- c(s, median_all_data.plot$samplenum[i], as.character(median_all_data.plot$time[i]), median_all_data.plot$category[i], (median_all_data.plot$nine_five_cluster[i]/unevolved)-1)
})

ninefive.fold.change.data <- data.frame(do.call(rbind, fold), stringsAsFactors = F)
colnames(ninefive.fold.change.data) <- append(colnames(median_all_data.plot[1:4]), "fold.data")
ninefive.fold.change.data$fold.data <- as.numeric(ninefive.fold.change.data$fold.data)

ninefive.fold.change.data$time <- factor(ninefive.fold.change.data$time, levels = c("t0", "t60", "t120", "t180"))
ninefive.fold.change.data$strain <- factor(ninefive.fold.change.data$strain, levels = c("GFP2", "GFPH", "myoH", "myoG")) ##

d5 <- ggplot(ninefive.fold.change.data, aes(time, fold.data)) + geom_boxplot(aes(col = strain), outlier.shape = NA) + 
  geom_jitter(width = 0.1) + geom_hline(yintercept = 0, linetype = "dotted") + facet_wrap(~strain, ncol = 2) + 
  theme(plot.title = element_text(size = 22), legend.position='right', strip.background = element_blank(),
        strip.text.x = element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=18), axis.title = element_text(size =18)) + 
  ylab("Fold change in cluster diameter") + xlab("Time") + scale_color_discrete(name="Strains") + ggtitle("95th-Percentile Cluster Diameter of Evolved Populations")

d5
### 5th percentile ###
fold <- lapply(1:nrow(median_all_data.plot), FUN = function(i) {
  s <- as.character(median_all_data.plot$strain[i])
  unevolved <- median_all_data.plot$five_cluster[which(median_all_data.plot$strain == s & median_all_data.plot$samplenum == 0)]
  fold.change <- c(s, median_all_data.plot$samplenum[i], as.character(median_all_data.plot$time[i]), median_all_data.plot$category[i], (median_all_data.plot$five_cluster[i]/unevolved)-1)
})

fifth.fold.change.data <- data.frame(do.call(rbind, fold), stringsAsFactors = F)
colnames(fifth.fold.change.data) <- append(colnames(median_all_data.plot[1:4]), "fold.data")
fifth.fold.change.data$fold.data <- as.numeric(fifth.fold.change.data$fold.data)

fifth.fold.change.data$time <- factor(fifth.fold.change.data$time, levels = c("t0", "t60", "t120", "t180"))
fifth.fold.change.data$strain <- factor(fifth.fold.change.data$strain, levels = c("GFP2", "GFPH", "myoH", "myoG")) ##

d6 <- ggplot(fifth.fold.change.data, aes(time, fold.data)) + geom_boxplot(aes(col = strain), outlier.shape = NA) + 
  geom_jitter(width = 0.1) + geom_hline(yintercept = 0, linetype = "dotted") + facet_wrap(~strain, ncol = 2) + 
  theme(plot.title = element_text(size = 22), legend.position='right', strip.background = element_blank(),
        strip.text.x = element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=18), axis.title = element_text(size =18)) + 
  ylab("Fold change in cluster diameter") + xlab("Time") + scale_color_discrete(name="Strains") + ggtitle("5th-Percentile Cluster Diameter of Evolved Populations")

d6

### combine percentile data ###
fifth.fold.change.data$percentile <- as.factor(rep("5th Percentile", length(nrow(fifth.fold.change.data))))
fifty.fold.change.data$percentile <- as.factor(rep("50th Percentile", length(nrow(fifty.fold.change.data))))
ninefive.fold.change.data$percentile <- as.factor(rep("95th Percentile", length(nrow(ninefive.fold.change.data))))
percentile.data <- data.frame(do.call(rbind, list(fifth.fold.change.data,fifty.fold.change.data,ninefive.fold.change.data)), stringsAsFactors = F)
mean.percentile.data <- percentile.data %>% group_by(strain, time, percentile) %>% summarise_at(vars(fold.data),list(fold.mean = mean)) %>% as.data.frame()

mean.percentile.data$percentile <- factor(mean.percentile.data$percentile, levels = c("5th Percentile", "50th Percentile", "95th Percentile"))

# plot
# single plot
#### adding to dataframe to plot as lines ###
df <- percentile.data
test <- percentile.data[which(percentile.data$time == "t0"),c(1,3:6)]
test1 <- do.call("rbind", replicate(5, test, simplify = FALSE)) %>% cbind(data.frame(samplenum = rep(c(1:5),each = 12)))
df <- rbind(df, test1)
df2 <- df[-c(which(df$time == "t0" & df$samplenum == 0)),]

leg <- c(expression(5^th ~ percentile), expression(50^th ~ percentile), expression(95^th ~ percentile))
x_labels <- c("0", "60", "120", "180")

df2$strain <- factor(df2$strain, levels = c("GFP2", "GFP2 + HEM3", "Myohemerythrin", "Myoglobin + HEM3"))

line.plot <- ggplot(df2, aes(time, fold.data, col= strain, lty = percentile, group = interaction(strain, percentile, samplenum))) + 
  geom_line(lwd=2) + scale_linetype_manual(values=c(3,1,6), labels = leg) + 
  geom_hline(yintercept=0 , linetype = "solid") + 
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        axis.text=element_text(size=34, color = "black"), legend.title = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1), 
        axis.text.x = element_text(size = 34, vjust = 0.5),
        axis.title.y = element_text(vjust=4), axis.title.x = element_text(vjust = -2), 
        legend.key = element_rect(fill = NA, color = NA), legend.position = "top",
        axis.ticks.y=element_line(size = 1), axis.ticks.x=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
        plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1,
        legend.box = "vertical", legend.justification = c(0.75,1)) + 
  ylab("Fold change in\n cluster diameter") + xlab("Days of evolution") + 
  scale_color_manual(values = c(denim,terra,leaf,purple)) + scale_x_discrete(labels = x_labels, expand = c(0,0)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=3))) 

ggsave(filename = "foldchange_v1.pdf", plot = line.plot, width = 12, height = 10, units = "in")


### anova and Tukey-hsd
df2 %>% group_by(percentile) %>% tukey_hsd(fold.data ~ strain)
aov.apply <- lapply(list.percentile, function(x){y <- x %>% aov(x$fold.data ~ x$strain, data=.) %>% summary()})
aov.results <- do.call("rbind", aov.apply) ## not significant (0.8, 0.24, 0.72)


b <- df2[which(df2$percentile == "5th Percentile"),]
summary(aov(fold.data ~ strain, data = b))
aov(fold.data ~ strain, data = b) %>% tukey_hsd()


### plotting with points
# single.plot <- ggplot(mean.percentile.data, aes(time, fold.mean, col=percentile)) + geom_point(size = 6) + 
#   geom_hline(yintercept=0 , linetype = "dashed") + 
#   theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
#         axis.text=element_text(size=34, color = "black"), legend.title = element_blank(),
#         plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1), 
#         axis.text.x = element_text(size = 34, vjust = 0.5),
#         axis.title.y = element_text(vjust=4), axis.title.x = element_text(vjust = -2), 
#         legend.key = element_rect(fill = NA, color = NA), legend.position = "top",
#         axis.ticks.y=element_line(size = 1), axis.ticks.length=unit(.25, "cm"),
#         plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1,
#         legend.justification = c(0.75,1)) + 
#   ylab("Average fold change\n in cluster diameter") + xlab("Days of evolution") + 
#   scale_color_manual(values = c(denim,terra,leaf), labels = leg) + scale_x_discrete(labels = x_labels) + 
#   annotate("point", x = 1, y = 0, colour = "black", size=9, shape = 18) # + labs(color = "Percentile")  
# 
# ggsave(filename = "mean_foldchange_v1.pdf", plot = single.plot, width = 12, height = 10, units = "in")

### four plots
#ggplot(mean.percentile.data, aes(time, fold.mean, col=percentile, group=strain)) + geom_point() + facet_wrap(vars(strain), nrow = 2)


#### two-way anova on cluster volume raw data
# all data
z2 <- all_data
z2$time <- as.numeric(z2$time)
z2$category <- factor(z2$category, levels = c("Control", "O2-Carriers"))
#summary(aov(volume_data ~ category * time, data = z2))
#TukeyHSD(aov(volume_data ~ strain * time, data = z2), which = "strain")
#z2 %>% anova_test(volume_data ~ category*time) # interaction between two independent variables
w. <- z2 %>% anova_test(volume_data ~ strain*time) %>% as.data.frame()
w <- aov(volume_data ~ strain*time, data = z2) %>% tukey_hsd() %>% as.data.frame()
w1 <- w[c(14,29,68,79,106,113,128,131),]

cluster.vol.aovtab <- ggtexttable(w., rows = NULL, theme = ttheme("light"))
cluster.vol.hsdtab <- ggtexttable(w1, rows = NULL, theme = ttheme("light"))

ggsave(filename = "cluster_vol_aov_tab.pdf", plot = cluster.vol.aovtab, width = 5, height = 4, units = "in")
ggsave(filename = "cluster_vol_hsd_tab.pdf", plot = cluster.vol.hsdtab, width = 10, height = 4, units = "in")

all_data_mean <- group_by(z2, strain, time) %>%
  summarise(
    count = n(),
    mean = mean(volume_data, na.rm = TRUE),
    sd = sd(volume_data, na.rm = TRUE)
  ) %>% as.data.frame()

o2.model <- lm(volume_data ~ strain + time + strain:time, data = z2)
anova(o2.model)

### ancova
aov(volume_data ~ strain + time + strain:time, data = z2) %>% tukey_hsd()

summary(aov(mean ~ strain + time + strain:time, data = all_data_mean))
summary(aov(pop.avg ~ strain + time + strain:time, data = z1)) # average of particles in each rep population


#### plotting raw data using kernal density ###
raw_data <- do.call("rbind", list.alldata)
ggplot(raw_data, aes(x = diameter_data, col = samplenum)) + geom_density() + facet_nested(time~category + strain)


d <- density(raw_data$diameter_data, bw = 0.5)

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

modes(d)


### finding peaks  by hand ###
evolved_peaks_file <- file.choose() %>% read.csv()
evolved_peaks_file$samplenum <- as.factor(evolved_peaks_file$samplenum)
evolved_peaks_file$radius <- ((evolved_peaks_file$volume)*(0.75/pi))^(1/3)

ggplot(evolved_peaks_file, aes(y=strain, x=radius, col = size)) + geom_point(size = 2)


### mixed model
# raw_data <- do.call("rbind", list.alldata)
# test1 <- raw_data[which(raw_data$strain == "Myoglobin + HEM3" & raw_data$samplenum == 2 & raw_data$time == 180),]
# ggplot(test1, aes(x=volume_data)) + geom_histogram(bins = 100, aes(y = stat(count / sum(count)))) + scale_x_log10()
# 
# plot_mix_comps <- function(x, mu, sigma, lam) {
#   lam * dnorm(x, mu, sigma)}
# 
# 
# 
# mixmdl1 <- normalmixEM(test1$diameter_data, k = 2)
# data.frame(x = mixmdl1$x) %>%
#   ggplot() +
#   geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
#                  fill = "white") +
#   stat_function(geom = "line", fun = plot_mix_comps,
#                 args = list(mixmdl1$mu[1], mixmdl1$sigma[1], lam = mixmdl1$lambda[1]),
#                 colour = "red", lwd = 1.5) +
#   stat_function(geom = "line", fun = plot_mix_comps,
#                 args = list(mixmdl1$mu[2], mixmdl1$sigma[2], lam = mixmdl1$lambda[2]),
#                 colour = "blue", lwd = 1.5) + xlim(c(0,100)) + geom_vline(xintercept = c(mixmdl1$mu[1], mixmdl1$mu[2]))



data.180.0 <- weighted_data.plot[which(weighted_data.plot$time == c(0,180)),] %>% group_split(samplenum, strain)

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# finding two peaks
y <- lapply(data.180.0, function(x){
  part <- x %>% filter(avg <= 5.5 & avg >= 3.3) %>% as.data.frame()
  peaks <- find_peaks(part$weight, m=8)
  left <- part[peaks[1],]
  right <- part[peaks[2],]
  #part1 <- part$weight
  #minima_peak <- find_peaks(-part$weight, m=50)
  #minima <- part[minima_peak[1],]
  rows <- rbind(left, right) %>% as.data.frame() #minima
  both <- cbind(rows, peak = c("left", "right")) # "minima"
})

y1 <- do.call("rbind", y) %>% as.data.frame()
y1 <- na.omit(y1)
y1$peak[1:4] <- rep("ancestor", 4)
gfph.5 <- y1[39,]

data.180 <- weighted_data.plot[which(weighted_data.plot$time == 180),] %>% group_split(samplenum, strain)

z <- lapply(data.180, function(x){
  test <- x %>% filter(avg <= 4.6 & avg >= 4) %>% as.data.frame()
  test1 <- test %>% summarise(test, diff = max(weight) - weight)
  #left.max <- test1[which.max(test1$weight),]
  df <- tail(test1[order(test1$diff, decreasing= F),], n = 10) %>% as.data.frame()
  m <- df[which.max(df$avg),]
  x1 <- x %>% as.data.frame()
  m1 <- match(m[1,1], x1[,1])
  left.area <- sum(x$weight[1:m1])
  right.area <- 1-left.area
  final <- cbind(m, m1, left.area, right.area)
  #left <- test[which.max(test$weight),]
  #r.peak <- x[nrow(test)+1:nrow(x),]
  #right <- test[which.max(r.peak$weight),]
  #final <- rbind(left, right)
  })
z1 <- do.call("rbind", z) %>% as.data.frame()
z1[,9] <- "minima"
z1 <- z1 %>% rename("peak" = "diff")
z1 <-z1[-18,]

y2 <- rbind(z1[,1:9], y1) #

#plotting peak points on top of cluster distributions
ggplot(anc.t180.df, aes(x = avg, y = weight, col = samplenum)) + geom_line(position = "identity", lwd = 1) + 
  facet_wrap(vars(strain), nrow = 2) + scale_x_continuous(labels = c(expression("10"^2),expression("10"^3),expression("10"^4),expression("10"^5),expression("10"^6))) + 
  #scale_x_continuous(label = scientific_10(breaks)) +
  #scale_x_continuous(labels  = function(x) format(breaks, scientific = T)) + 
  ylab("Frequency") + xlab(expression(paste("Cluster Volume (", mu,"m"^{3},")"))) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(strip.background = element_blank(), legend.key = element_blank(),
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=34), 
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), legend.position = "none", 
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_color_manual(name = "Replicate Populations", values = c("black", terra, denim, leaf, beige, purple)) + geom_point(data = y1, aes(x = avg, y = weight), size = 3) + 
  geom_vline(data = z1, aes(xintercept = avg, col = samplenum), linetype = "dotted")



z1.ordered <- z1[order(z1$samplenum, decreasing = F),]
z1.area <- z1.ordered[,11:12]
a <- as.vector(t(z1.area)) %>% as.data.frame()
colnames(a) <- "area"
empty <- data.frame("area" = rep(1, 4))
a <- rbind(empty, a)
y1 <- y1[-39,]
y1$area <- a
gfph.5$area <- as.numeric(1)
y1 <- rbind(y1, gfph.5)

#y2 <- y1[5:nrow(y1),]
y1$area <- as.numeric(unlist(y1$area))

y2 <- y1[which(y1$category == "O2-Carriers"),]

# plotting volume at peak versus population proportion, this can somewhat show spread
ggplot(y1, aes(x = avg, y = area, col = samplenum, group = strain)) + geom_point(size = 5) + geom_line(aes(group = samplenum), lwd = 1.5) + 
  facet_wrap(vars(strain), nrow = 1) + scale_x_continuous(labels = c(expression("10"^3.5),expression("10"^4.0),expression("10"^4.5),expression("10"^5))) + 
  scale_color_manual(name = "Replicate Populations", labels = c('Ancestor', '1', '2', '3','4','5'), values = c("black", terra, denim, leaf, beige, purple)) + 
  xlab(expression(paste("Cluster Volume (", mu,"m"^{3},")"))) + ylab("Proportion of Population") 

# plotting fold difference in cluster volume peaks versus difference in population proportion
volume.area <- y1[5:42,c(3,5:8,10)]
# right peak/left peak, right peak proportion - left peak proportion
t <- volume.area %>% group_by(strain, samplenum) %>% summarise(volume = last(avg)/first(avg), area = last(area)/first(area)) %>% as.data.frame()

ggplot(t, aes(x = volume, y = area, col = samplenum, group = strain)) + geom_point(size = 5) + 
  facet_wrap(vars(strain), nrow = 2) +
  scale_color_manual(name = "Replicate Populations", values = c(terra, denim, leaf, beige, purple)) + 
  xlab("Fold difference in cluster volume") + ylab("Fold difference in proportion of population") + 
  ggtitle("Fold differences in volume and population proportion between larges and smalls in a population")

# plotting fold difference of each peak's volume size relative to ancestor
s <- y1 %>% group_by(strain) %>% summarise(samplenum, volume = avg/first(avg), area) %>% as.data.frame() %>% subset(., samplenum!=0)

ggplot(s, aes(x = volume, y = area, col = samplenum, group = strain)) + geom_point(size = 5) + geom_line(aes(group = samplenum), lwd = 1.5) + 
  facet_wrap(vars(strain), nrow = 2) +
  scale_color_manual(name = "Replicate Populations", values = c(terra, denim, leaf, beige, purple)) + 
  xlab("Fold difference in cluster size relative to ancestor") + ylab("Proportion of Population") +
  geom_vline(xintercept = 1, linetype = "dotted")

# x-axis of the above plot

s1 <- y1 %>% group_by(strain) %>% summarise(samplenum, category, peak, volume = avg/first(avg)) %>% as.data.frame() %>% subset(., samplenum!=0)
s1$strain <- factor(s1$strain, levels = c("GFP2", "Myohemerythrin", "GFP2 + HEM3", "Myoglobin + HEM3"))

ggplot(s1, aes(x= volume, y= strain)) +
  geom_line(aes(group = strain))+
  geom_point(aes(color=samplenum), size=7, shape=21, stroke=3, position=position_jitter(h=0.08)) + 
  geom_vline(xintercept = 1, linetype = "dotted") + scale_color_manual(name = "Populations", values = c(terra, denim, leaf, beige, purple)) +
  xlim(c(0.82,1.26)) + xlab("Fold difference in cluster \nsize relative to ancestor") + 
  theme(axis.title = element_text(size =34), axis.title.y = element_blank(), legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=15), legend.title = element_text(size=18), axis.text=element_text(size=34, color = "black"), 
        plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), axis.title.x = element_text(vjust = -0.5),
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5))

# finding average of datapoints left and right of 1.0
s2 <- rbind(s1, cbind(s1[19, 1:3], peak = "right", volume = s1[19,5]))
s3 <- s2 %>% group_by(strain, peak) %>% summarise(average = mean(volume)) %>% as.data.frame()


ggplot(s1, aes(x= volume, y= strain)) +
  geom_line(aes(group = strain))+
  geom_point(aes(color=samplenum), size=7, shape=21, stroke=3, position=position_jitter(h=0.08)) + 
  geom_vline(xintercept = 1, linetype = "dotted") + scale_color_manual(name = "Populations", values = c(terra, denim, leaf, beige, purple)) +
  xlim(c(0.82,1.26)) + xlab("Fold difference in cluster \nsize relative to ancestor") + 
  theme(axis.title = element_text(size =34), axis.title.y = element_blank(), legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=15), legend.title = element_text(size=18), axis.text=element_text(size=34, color = "black"), 
        plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("white"), axis.title.x = element_text(vjust = -0.5),
        axis.text.x = element_text(size = 30, vjust = 0.5, hjust = 0.5)) + geom_point(data = s3, aes(x = average, y = strain, shape = peak), size = 5, show.legend = F)

s4 <- cbind(s2[order(s2$strain, decreasing = F),], group = c(rep("one", 20), rep("two", 20)))
S <- s4 %>% group_by(peak, group) %>% t_test(data =., volume ~ strain) %>% add_significance("p") %>% add_y_position() %>% as.data.frame()


# Separating peak data for Will
w180 <- lapply(data.180, function(x){
  test <- x %>% filter(avg <= 4.6 & avg >= 4) %>% as.data.frame()
  test1 <- test %>% summarise(test, diff = max(weight) - weight)
  df <- tail(test1[order(test1$diff, decreasing= F),], n = 10) %>% as.data.frame()
  m <- df[which.max(df$avg),]
  x1 <- x %>% as.data.frame()
  m1 <- match(m[1,1], x1[,1])
  left.area <- sum(x$weight[1:m1])
  right.area <- 1-left.area
  final <- cbind(m, m1, left.area, right.area)
  
  weighted.180.f <- x %>% mutate(group.size = case_when(x[,3] < final[,3] ~ "Small", TRUE ~ "Large"))
  
})

w180.f <- do.call("rbind", w180)
w180.f <- w180.f[rep(row.names(w180.f), w180.f$n), c(1,3,5:9)]
w180.f$volume.size <- 10^(w180.f$avg)
w180.f[which(w180.f$strain == "GFP2 + HEM3" & w180.f$samplenum == 5),7] <- NA

write.csv(w180.f, "O2-evolved_groupsize-data.csv", row.names = F)

