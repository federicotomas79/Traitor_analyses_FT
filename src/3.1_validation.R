library(tidyverse)
library(DescTools)
library(ggplot2)
library(gridExtra)

#setwd("write here the path to your working directory") # setting the working directory folder

df = read.csv("DiasMorph_quantitative_traits.csv") # available at XXXXX

#names(df)

# reducing size of df. If working with a different version of DiasMorph, make sure it has identifying columns and length and width measurements
df = df[,c(1:10)] 

# overview of dataset
##samples per taxon
sp_hist = as.data.frame(table(df$scientificName))
hist(sp_hist$Freq, breaks=50)
length(sp_hist$Freq[sp_hist$Freq >170])
sp_hist$Freq2 <- ifelse(sp_hist$Freq > 170, 173, sp_hist$Freq) #for histogram plot only

#pdf(file = "histogram.pdf", width = 6.2, height = 4)
ggplot(sp_hist, aes(x=Freq2)) +
  geom_histogram(binwidth=5, color="white", alpha=0.9, boundary=0) +
  scale_x_continuous(breaks = seq(0, 175, by = 10), labels=c(seq(0,160, by=10), "170+")) +
  xlab("Number of sampled diaspores") +
  ylab("Number of taxa") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(hjust=0))
#dev.off()

# validation -------------------------------------------------------------------

# disregarding species with less than 5 reps for size comparison -------------
sp = table(df$image_name)
notEnough = names(sp)[sp < 5]
df = df[! (df$image_name %in% notEnough), ]

# removing species with more than one measurement ------------------------------
# (e.g., with and without specific appendages)
# increases the chances of comparing the same structures
colnames(df)
id_df = unique(df[,c(3:8)])
df_n_count = table(id_df$scientificName)
two_states = names(df_n_count)[df_n_count > 1]
df <- df[! (df$scientificName %in% two_states), ]

length(df$scientificName)

# summary of DiasMorph (Traitor) measurements ----------------------------------
names(df)

traitor_sum <- df %>%
  group_by(scientificName, scientificNameAuthorship, spec.name.ORIG, genus, family, missing_structures) %>%
  summarise(traitor_length = mean(length), traitor_length_sd = sd(length), traitor_width = mean(width), traitor_width_sd = sd(width))

print(traitor_sum)

# manual measurements ----------------------------------------------------------

manual = read.csv("manual_measurements.csv") # see paper's supplementary material
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

g1 = 
ggplot(data = data) + 
  geom_segment(aes(x=-1.3, y=10, xend=10, yend=10), linewidth=0.8, colour="grey70", linetype=3) +
  geom_segment(aes(y=-1, x=10, xend=10, yend=10), linewidth=0.8, colour="grey70", linetype=3) +
  geom_point(aes(x=manual_length , y = traitor_length), colour = "#0277BD", size = 1.2, alpha = 0.7) +
#  geom_errorbar(aes(x = manual_length, ymin = traitor_length - traitor_length_sd, ymax = traitor_length + traitor_length_sd), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for Traitor's measurements
#geom_errorbar(aes(y = traitor_length, xmin = manual_length - sd_manual_length, xmax = manual_length + sd_manual_length), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for manual measurements
  geom_abline(intercept = 0, slope = 1, colour = "black", linewidth=0.5, linetype=5) +
  scale_y_continuous(limits=c(-1,31.5),expand = c(0, 0)) +
  scale_x_continuous(limits=c(-1.3,33.3),expand = c(0, 0)) +
  ylab(paste("length Traitor (mm)")) +
  xlab(paste("length manual (mm)")) +
  annotate(geom = "text", label = "a", x = 1, y = 29, size = 5) +
  annotate(geom = "text", label = "rho==0.979", parse = TRUE, x = 23, y = 7, size = 5, hjust = 0) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
g1

g2 = 
  ggplot(data = data) + 
  geom_segment(aes(x=-1, y=5, xend=5, yend=5), linewidth=0.8, colour="grey70", linetype=3) +
  geom_segment(aes(y=-0.9, x=5, xend=5, yend=5), linewidth=0.8, colour="grey70", linetype=3) +
  geom_point(aes(x=manual_width , y = traitor_width), colour = "#0277BD", size = 1.2, alpha = 0.7) + 
#  geom_errorbar(aes(x = manual_width, ymin = traitor_width - traitor_width_sd, ymax = traitor_width + traitor_width_sd), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for Traitor's measurements
#geom_errorbar(aes(y = traitor_width, xmin = manual_width - sd_manual_width, xmax = manual_width + sd_manual_width), width = 0, alpha = 0.7, colour = "#0277BD", size =0.3) + # add error bar with sd for manual measurements
  geom_abline(intercept = 0, slope = 1, colour = "black", size =0.5, linetype=5) +
  scale_y_continuous(limits=c(-0.9,23.7),expand = c(0, 0)) +
  scale_x_continuous(limits=c(-1.1,24.7),expand = c(0, 0)) +
  ylab(paste("width Traitor (mm)")) +
  xlab(paste("width manual (mm)")) +
  annotate(geom = "text", label = "c", x = 0.65, y = 21.5, size = 5) +
  annotate(geom = "text", label = "rho==0.984", parse = TRUE, x = 17, y = 5, size = 5, hjust = 0) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
g2

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

# genus and family ------------------------------
names(data)
taxa_fin = unique(data[,c(1,4,5)])
length(unique(taxa_fin$genus))
length(unique(taxa_fin$family))

taxa_fin %>% count(family, sort = TRUE)
