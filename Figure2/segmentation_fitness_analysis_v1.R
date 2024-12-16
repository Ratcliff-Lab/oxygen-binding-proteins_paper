### Sorting particle data from ImageJ/Fiji, then calculates relative fitness ###

###history###
# from 'comp5.R'# script
# Version: v6 - renamed

### what does this script do?###
#1. Allows user to interactively select folder to read and analyze particle segmented ImageJ .csv outputs
  #1a. Reads .csv, sorts the 'Mean' column, filters out 10% of data, creates indices to find cutoff of small and large values to differentitate GFP vs Non-GFP, calculates biomass, calculates biomass, and proportion of two populations
#2. Output from function is a list, reformatted into a dataframe with new column names, adds metadata of experimental groups
#3. Finds relative fitness
#4. Generates box plot and a .csv file for plotting dot plot in another script

### required directory organization ###
# None

### required filenaming and columnnaming scheme ###
# csv filenames:  '#3_a_t0_1.csv', '#3_e_t0_2.csv', '#4_a_t0_1.csv', '#4_e_t3_2.csv'
# csv ccolumn contents:  'Area'	'Mean' 'Min' ;'Max'

### referenece links for help ###
# https://lenski.mmg.msu.edu/ecoli/srvsrf.html

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

# choosing directory 
wd <- choose.dir(default = NA, caption = "Select Folder")
setwd(wd)

f_day <- 3

# reads all raw ImageJ files with .csv extension in directory, put only these files in folder
filenames_bf <- list.files(path = wd, pattern = ".csv", full.names=TRUE, recursive = FALSE)


# long (probably not necessary) function to sort, apply filters, and apply cutoff to differentitate fluorescent/nonfluorescent particles
# assigns output to process_intensity object; output is a list
process_intensity <- lapply(1:length(filenames_bf), FUN = function(i) {
  # extract file names 
  file.name <- tools::file_path_sans_ext(basename(filenames_bf[[i]]))
  x <- str_split(file.name, pattern = "_")
  x1 <- x[[1]][c(1,2,3)] #c(1,2,3)
  file.name.without.techrep <- str_c(x1,collapse='_') 
  
  results.data <- read.csv(filenames_bf[[i]]) %>% data.frame # read .csv files
  sorted.results.data <- results.data[order(results.data$Mean), c(2,3)] # increasingly sort Mean column 
  rownames(sorted.results.data) <- c(1:nrow(sorted.results.data)) # rename row names 
  diff.sorted.results.data <- data.frame(Mean.Intensities = diff(sorted.results.data$Mean)) # differences by 1
  diff.sorted.results.data$x <- 1:nrow(diff.sorted.results.data) # a column of indices 
  ###p <- ggplot(diff.sorted.results.data, aes(x = x, y = Mean.Intensities)) + geom_line() # plot of differences
  
  cut.upper.tenth <- round(0.9*nrow(sorted.results.data)) # getting indices of first 90% of data
  sec.diff.sorted.results.data <- data.frame(Diff.Mean.Intensities = diff(diff.sorted.results.data$Mean.Intensities[1:cut.upper.tenth])*-1) # second round of differences of 90% of data
  sec.diff.sorted.results.data$x <- 1:nrow(sec.diff.sorted.results.data) # a column of indices 
  cutoff.index <- match(max(sec.diff.sorted.results.data$Diff.Mean.Intensities), sec.diff.sorted.results.data$Diff.Mean.Intensities) # find index of cutoff
  
  total.biomass <- sum((sorted.results.data$Area)**(3/2)) # total biomass
  Non.GFP <- sum((sorted.results.data$Area[1:cutoff.index])**(3/2)) # total biomass of non-GFP-tagged strain
  Non.GFP.prop <- sum((sorted.results.data$Area[1:cutoff.index])**(3/2))/total.biomass # proportion of non-GFP-tagged strain
  avg.Non.GFP <- mean((sorted.results.data$Area[1:cutoff.index])**(3/2)) # average biomass of non-GFP-tagged strain
  
  c(file.name, file.name.without.techrep, Non.GFP, total.biomass, Non.GFP.prop, avg.Non.GFP)
}) 

# bind elements in list into a readable dataframe
process_intensity <- as.data.frame(do.call(rbind, process_intensity), stringsAsFactors = F)
# renaming column names
colnames(process_intensity) <- c("File.Name", "File.Name1", "Nonfluorescent.Biomass", "Total.Biomass","Nonfluorescent.Prop", "Average.Area")
# converting characters into numerics
process_intensity$Nonfluorescent.Biomass <- as.numeric(process_intensity$Nonfluorescent.Biomass)
process_intensity$Total.Biomass <- as.numeric(process_intensity$Total.Biomass)
process_intensity$Average.Area<- as.numeric(process_intensity$Average.Area)

# summing technical replicates
processed.intensity <- aggregate(process_intensity[,3:4], by = list(File.Name1 = process_intensity$File.Name1), FUN = sum)
# normalizing nonfluorescent biomass to total population biomass
processed.intensity$Nonfluorescent.Prop <- processed.intensity$Nonfluorescent.Biomass/processed.intensity$Total.Biomass

### processed.intensity <- process_intensity ###

# adding columns of metadata
for(i in 1:nrow(processed.intensity)) {
  #processing data
  comp_group <- str_extract(processed.intensity$File.Name1[i], pattern = '(?<=#)(\\d)+') 
  processed.intensity$Bio.Replicate[i] <- str_extract(processed.intensity$File.Name[i], pattern = '(?<=_)[a-zA-Z]') #%>% as.integer
  processed.intensity$Day[i] <- str_extract(processed.intensity$File.Name[i], pattern = '(?<=t)(\\d)+') %>% as.integer
  processed.intensity$Competing.Groups[i] <- switch(comp_group %>% as.integer(), "GFP2 - GOB8","GFP2+HEM - GOB8", "GFP2 - myoH", "GFP2+HEM - myoG+HEM") #"space"
  processed.intensity$Variable[i] <- switch(comp_group, "1" = "Control", "2" ="Control", "O2-Carrier")
  processed.intensity$Phenotype[i] <- switch(comp_group %>% as.integer(), "Normal","Normal", "Normal", "Normal", "Normal")
}

# reorganizing factors
processed.intensity$Competing.Groups <- factor(processed.intensity$Competing.Groups, levels = 
                                                 c("GFP2 - GOB8", "GFP2+HEM - GOB8", "GFP2 - myoH", "GFP2+HEM - myoG+HEM", "space"))


############# relative fitness ##################

# grabbing specific columns based on position number
rel.data <- processed.intensity[which(processed.intensity$Day == 0 | processed.intensity$Day == f_day), c(1, 4:9)] #day
#rel.data$Competing.Groups <- factor(levels(rel.data$Competing.Groups)[1:4])
rel.data$Bio.Replicate <- as.factor(rel.data$Bio.Replicate)

levels.c.g <- length(levels(rel.data$Competing.Groups)[1:4])
levels.b.r <- length(levels(rel.data$Bio.Replicate))
total.levels <- levels.c.g*levels.b.r

prop.col <- 2 # 2

relative.data <- list()
subdata <- list()
for (i in 1:levels.c.g) {
  comp.group <- levels(rel.data$Competing.Groups)[i]
  for (j in 1:levels.b.r) {
    bio.rep <- levels(rel.data$Bio.Replicate)[j]
    initial.day <- as.numeric(rel.data[which(rel.data$Competing.Groups == comp.group & rel.data$Day == 0 & rel.data$Bio.Replicate == bio.rep), prop.col])
    final.day <- as.numeric(rel.data[which(rel.data$Competing.Groups == comp.group & rel.data$Day == f_day & rel.data$Bio.Replicate == bio.rep), prop.col]) #day
    #fitness <- log((final.day - (final.day*initial.day))/(initial.day - (final.day*initial.day)))
    fitness <- log((final.day/initial.day)*(200**f_day))/log((((1-final.day)/(1-initial.day))*(200**f_day))) - 1 #day
    subdata[[j]] <- c(comp.group, bio.rep, fitness)
  }
  relative.data[[i]] <- subdata
}

fitness.data <- as.data.frame(matrix(unlist(relative.data[3:4]), ncol=3, byrow = TRUE)) # looking at only oxygen-carriers
colnames(fitness.data) <- c("Competing.Group", "Bio.Replicate", "Relative.Fitness")

for (i in 1:nrow(fitness.data)) {
  comp_group <- sub(".*- ", "", fitness.data$Competing.Group[i])
  fitness.data$Variable[i] <- switch(comp_group, "GOB8" = "Control", "myoH" ="O2-Carrier", "O2-Carrier")
}

# final datapoints for plotting
fitness.data$Relative.Fitness <- as.numeric(as.character(fitness.data$Relative.Fitness))

### when analyzing different oxygen conditions, changing the following column ###
fitness.data$Treatment <- "Low O2"

# median for boxplot
mid <- aggregate(Relative.Fitness ~  Competing.Group, fitness.data, median)
mid$Relative.Fitness <- round(mid$Relative.Fitness, 4)

# colors for plotting
blue <- "#B7E2F7"
green <- "#92E0D9"
gold <- "#FFC654"
lavender <- "#EBE8FF"
pink <- "#FF545A"

r.f. <- ggplot(fitness.data, aes(x = Competing.Group, y = Relative.Fitness)) + geom_point(aes(col = Competing.Group), size = 0, shape = 16, show.legend = T) + geom_boxplot(aes(col = Competing.Group), show.legend = F) + geom_jitter(width = 0.1) + #aes(col = Competing.Group),
  geom_hline(yintercept=0, linetype = "dashed") + theme(strip.background = element_blank(),
                                                        strip.placement = "outside", strip.text = element_text(size=15),
                                                        legend.key.size = unit(1, 'cm'), plot.title = element_text(size = 22),
                                                        axis.title = element_text(size =18), legend.text = element_text(size=12),
                                                        legend.title = element_text(size=15), axis.text=element_text(size=12),
                                                        plot.margin = margin(20, 150, 20, 150), panel.background = element_rect("white"),
                                                        axis.line = element_line(color = "black"), panel.grid.major = element_line(color = 'gray',
                                                                                                                                   linetype = 'dotted'), axis.text.x=element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values = c(gold, pink)) + ylab(expression(paste('Relative fitness (W) of '~O[2]~'-Carrier Strain'))) + labs(col = "Competing Groups", caption = "2021-6-9") + ggtitle("Relative Darwinian Fitness of "~O[2]~"-Carriers") +
  theme(plot.caption.position = "panel", plot.tag.position = "bottom", plot.caption = element_text(hjust = 0, vjust = -2, size = 20, color = "gray"), legend.key = element_rect(fill = NA)) + scale_y_continuous(limits=c(-0.01,0.04)) + 
  geom_text(data = mid, aes(label = Relative.Fitness, y = Relative.Fitness, col = Competing.Group), vjust = 2, size = 5, show.legend = F) + guides(color = guide_legend(override.aes = list(size = 8)))


r.f.

# exporting values for plotting dot-plot in another script
write.csv(fitness.data, paste(str_split(wd, "/")[[1]][7], "_fitness.data_20241006", ".csv", sep =''), row.names = F)


