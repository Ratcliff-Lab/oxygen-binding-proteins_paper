### plotting selective advantage model results ###

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
library(RColorBrewer)
library(tidyr)
library(pals)

# choosing files
wd <- choose.dir(default = NA, caption = "Select Folder")
sel.adv <- choose.dir(default = NA, caption = "Select Folder")
setwd(wd)

sel.adv.data <- file.choose() %>% read.csv() 
  # "/Tonys_modeldata/selective_advantage/Myo70COMPLETEwithlowO2growthadvantage.csv"
sel.adv.data$Radius <- as.factor(sel.adv.data$Radius)
sel.adv.data$Myoglobin <- as.factor(sel.adv.data$Myoglobin)

low.sel.adv.data <- file.choose() %>% read.csv() 
  # "/Tonys_modeldata/selective_advantage/lowoxyfillin.csv"
low.sel.adv.data$Radius <- as.factor(low.sel.adv.data$Radius)
low.sel.adv.data$Myoglobin <- as.factor(low.sel.adv.data$Myoglobin)

# getting the "advantage" column
low.sel.adv.data <- low.sel.adv.data %>% group_by(Radius) %>%mutate(Advantage = (Metabolism/Metabolism[Myoglobin == 0])-1)

sel.data <- rbind(sel.adv.data, low.sel.adv.data)

filter.sel.adv.data <- 
  sel.data[which(sel.data$Radius %in% c(5, 10, 20) & sel.data$Myoglobin %in% c(0.01, 0.05, 0.1, 0.2)),]

purple <- "#7678ed"
leaf <- "#81B29A"
beige <- "#F2CC8F"
terra <- "#E07A5F"
denim <- "#457b9d"

filter.sel.adv.data$Radius <- factor(filter.sel.adv.data$Radius, labels = c(bquote(5~mu*m), bquote(10~mu*m), bquote(20~mu*m)))

s.a.l <- ggplot(filter.sel.adv.data, aes(x = Oxygen, y = Advantage, group = interaction(Myoglobin, Radius),
                         col = Myoglobin)) + geom_line(lwd = 1.5) +
  geom_smooth(se=FALSE, col = NA) + scale_x_continuous(breaks = c(0.05, 0.15, 0.25), expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(0, 0.55, by = 0.1), expand = c(0, 0)) + 
  facet_wrap2(~Radius, labeller = label_parsed, axes = "all", remove_labels = "all") +
  theme(strip.background = element_blank(), 
        strip.placement = "outside", strip.text = element_text(size=34, vjust = 1),
        legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.65,1), legend.title = element_text(size=34), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1.3,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top", legend.background = 
          element_rect(fill="transparent", size=0.5),
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5, angle = 0), panel.spacing.x = unit(2.5, "lines"), 
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio = 1) +
  scale_linetype_manual(labels = c(expression(paste("5 ",mu,"m")), expression(paste("10 ",mu,"m")),expression(paste("20 ",mu,"m"))), 
                        values = c('dotted','dashed','solid')) + 
  scale_color_manual(name = "Myoglobin", labels = c("0.01 mM", "0.05 mM", "0.1 mM", "0.2 mM"), values = c(leaf, beige, terra, purple)) + 
  xlab("Oxygen (mM)") + ylab("Selective Advantage")

ggsave(filename = "selective_adv_lines_v8.pdf", plot = s.a.l, width = 15, height = 10, units = "in")


### heatmap
sel.adv.data.h <- file.choose() %>% read.csv(., row.names = 1, header = T, check.names = F) # arguments to include col and row names
  # "/Tonys_modeldata/selective_advantage/HighResHeatMap_trans.csv"
sel.adv.data.h <- sel.adv.data.h[,1:14]
sel.adv.data.h <- tibble::rownames_to_column(sel.adv.data.h, "Oxygen")

long <- sel.adv.data.h %>% 
  pivot_longer(
    cols = `5`:`70`, 
    names_to = "Radius",
    values_to = "Globin.effects"
  ) %>% mutate_at(., c("Radius"), as.numeric)

###s.a.d <- long[which(long$Oxygen > 0.009),]
#colors <- brewer.pal(9, 'YlOrRd')
#red <- "#ed1c1c"
#blue <- "#3fa4e8"

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')

# s.a.h <- ggplot(s.a.d, aes(Oxygen, Radius, fill= Globin.effects)) + 
#   geom_raster() + 
#   theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
#         legend.box = "vertical", legend.justification = c(0.5,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
#         axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
#         panel.background = element_rect("transparent"), legend.position = "top",
#         axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"), 
#         panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1.45) + 
#   scale_y_continuous(breaks = seq(0, 70, by = 10), expand = c(0, 0)) + scale_x_discrete(breaks = seq(0.05, 0.25, by = 0.1)) + 
#   scale_fill_gradient2(name = c("Selective\nAdvantage"), midpoint = 0.2, low = blue, high = red) + 
#   ylab(expression(paste("Radius ","(",mu, "m)"))) + xlab("Oxygen (mM)") #+ coord_fixed(ratio = 0.3) #0.3

# s.a.h <- ggplot(long, aes(Oxygen, Radius, fill= Globin.effects)) + 
#   geom_raster() + 
#   theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
#         legend.box = "vertical", legend.justification = c(0.5,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
#         axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
#         panel.background = element_rect("transparent"), legend.position = "top",
#         axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"),
#         panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1.45) + 
#   scale_y_continuous(breaks = seq(0, 70, by = 10), expand = c(0, 0)) + scale_x_discrete(breaks = seq(0.0, 0.25, by = 0.05)) + 
#   scale_fill_gradient(name = c("Selective\nAdvantage"), low = blue, high = red, guide = guide_colorbar(ticks.colour = "black")) + 
#   ylab(expression(paste("Radius ","(",mu, "m)"))) + xlab("Oxygen (mM)")

s.a.h <- ggplot(long, aes(Oxygen, Radius, fill= Globin.effects)) + 
  geom_raster() + 
  theme(legend.key.size = unit(1.5, 'cm'), axis.title = element_text(size =34), legend.text = element_text(size=30, margin = margin(r = 10, unit = "pt")), 
        legend.box = "vertical", legend.justification = c(0.5,1), legend.title = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), axis.title.y = element_text(vjust = 4), axis.title.x = element_text(vjust = -0.5),
        axis.text=element_text(size=34, color = "black"), plot.margin = unit(c(1,1,1,1), "cm"), axis.line = element_line(color = "black", size = 1),
        panel.background = element_rect("transparent"), legend.position = "top",
        axis.text.x = element_text(size = 34, vjust = 0.5, hjust = 0.5), panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(2, "lines"), plot.background = element_rect(fill = "transparent",colour = NA), aspect.ratio=1.45) + 
  scale_y_continuous(breaks = seq(0, 70, by = 10), expand = c(0, 0)) + scale_x_discrete(breaks = seq(0.0, 0.25, by = 0.05)) + 
  scale_fill_gradientn(name = c("Selective\nAdvantage"), colours = hm.palette(100), guide = guide_colorbar(ticks.colour = "black")) + 
  ylab(expression(paste("Radius ","(",mu, "m)"))) + xlab("Oxygen (mM)")

# scale_fill_gradient2(midpoint = 0.2, low = scales::muted(blue), high = scales::muted(terra), breaks = seq(0,0.4, by = 0.1))
# scale_fill_gradientn(colours = colors), scale_fill_distiller(palette = "RdBu") 
ggsave(filename = "selective_adv_heat_v11.pdf", plot = s.a.h, width = 15, height = 14, units = "in")


