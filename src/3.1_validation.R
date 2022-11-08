library(tidyverse)
library(DescTools)
library(ggplot2)
library(gridExtra)

#setwd("XXX") # insert working directory

df = read.csv("DiasMorph_v1.csv") # available at XXXXX

#names(df)

# reducing size of df. If working with a different version of DiasMorph, make sure it has identifying columns and length and width measurements
df = df[,c(1:10)] 

# validation -------------------------------------------------------------------

# disregarding species with less than 5 reps for size comparison -------------
sp = table(df$group)
notEnough = names(sp)[sp < 5]
df = df[! (df$group %in% notEnough), ]

# removing species with more than one measurement ------------------------------
# (e.g., with and without specific appendages)
# increases the chances of comparing the same structures
colnames(df)
id_df = unique(df[,c(3:8)])
df_n_count = table(id_df$scientificName)
two_states = names(df_n_count)[df_n_count > 1]
df <- df[! (df$scientificName %in% two_states), ]

# summary of DiasMorph (Traitor) measurements ----------------------------------
names(df)

traitor_sum <- df %>%
  group_by(scientificName, scientificNameAuthorship, spec.name.ORIG, genus, family, missing_structures) %>%
  summarise(traitor_length = mean(length), traitor_length_sd = sd(length), traitor_width = mean(width), traitor_width_sd = sd(width))

print(traitor_sum)

# manual measurements ----------------------------------------------------------

manual = read.csv("manual_measurements.csv") # available on data folder on GitHub and as supplementary material

names(manual)

# merge manual and Traitor datasets --------------------------------------------

data = merge.data.frame(traitor_sum, manual, by = c("scientificName", "scientificNameAuthorship", "spec.name.ORIG", "genus", "family"), all.x = FALSE, all.y = FALSE)

data = droplevels(data)
length(levels(as.factor(data$scientificName)))

# comparison Traitor vs. manual ------------------------------------------------
# Calculate Lin's concordance correlation coefficient

# length
ccc_len = CCC(data$manual_length, data$traitor_length, ci = "z-transform", conf.level = 0.95)
ccc_len$rho.c

# width
ccc_wid = CCC(data$manual_width, data$traitor_width, ci = "z-transform", conf.level = 0.95)
ccc_wid$rho.c

# plot results -----------------------------------------------------------------
# Validation figure

g1 = 
ggplot(data = data) + 
  geom_point(aes(x=manual_length , y = traitor_length), colour = "#0277BD", size = 1.2, alpha = 0.7) +
#  geom_errorbar(aes(x = manual_length, ymin = traitor_length - traitor_length_sd, ymax = traitor_length + traitor_length_sd), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for Traitor's measurements
#geom_errorbar(aes(y = traitor_length, xmin = manual_length - sd_manual_length, xmax = manual_length + sd_manual_length), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for manual measurements
 geom_abline(intercept = 0, slope = 1, colour = "black", size =0.5, linetype=5) +
  geom_vline(xintercept = 10, lwd=0.8, colour="grey80", linetype=3) +
  ylab(paste("length Traitor (mm)")) +
  xlab(paste("length manual (mm)")) +
  annotate(geom = "text", label = "a", x = 1, y = 29, size = 5) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
#g1

g2 = 
  ggplot(data = data) + 
  geom_point(aes(x=manual_width , y = traitor_width), colour = "#0277BD", size = 1.2, alpha = 0.7) + 
#  geom_errorbar(aes(x = manual_width, ymin = traitor_width - traitor_width_sd, ymax = traitor_width + traitor_width_sd), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for Traitor's measurements
#geom_errorbar(aes(y = traitor_width, xmin = manual_width - sd_manual_width, xmax = manual_width + sd_manual_width), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for manual measurements
  geom_abline(intercept = 0, slope = 1, colour = "black", size =0.5, linetype=5) +
  geom_vline(xintercept = 5, lwd=0.8, colour="grey80", linetype=3) +
  ylab(paste("width Traitor (mm)")) +
  xlab(paste("width manual (mm)")) +
  annotate(geom = "text", label = "c", x = 0.7, y = 21, size = 5) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
#g2

g3 = 
ggplot(data = data) + 
  geom_point(aes(x=manual_length , y = traitor_length), colour = "#0277BD", size = 1.2, alpha = 0.7) + 
  geom_abline(intercept = 0, slope = 1, colour = "black", size =0.5, linetype=5) +
  xlab(paste("length manual (mm)")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  annotate(geom = "text", label = "b (detailed view)", x = 2.5, y = 9.5, size = 5) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.title.y = element_blank())
#g3

g4 = 
ggplot(data = data) + 
  geom_point(aes(x=manual_width , y = traitor_width), colour = "#0277BD", size = 1.2, alpha = 0.7) + 
  geom_abline(intercept = 0, slope = 1, colour = "grey30", linetype=5) +
  xlab(paste("width manual (mm)")) +
  ylim(0, 5) +
  xlim(0, 5) +
  annotate(geom = "text", label = "d (detailed view)", x = 1.3, y = 4.8, size = 5) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.title.y = element_blank())
#g4

pdf(file = "validation.pdf", width = 8, height = 6)
grid.arrange(arrangeGrob(g1,g3,g2,g4, as.table= TRUE, ncol=2, nrow=2, heights=c(5,5), widths=c(6.5,4.5)))
dev.off()

