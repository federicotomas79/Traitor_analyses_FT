library(ggplot2)

#setwd("write here the path to your working directory") # setting the working directory folder

df = read.csv("DiasMorph_v1.csv")  # available at XXXXX
# head(df)

names(df)

# reducing size of df. If working with a different version of DiasMorph, make sure it has identifying columns and colour measurements
df = df[,c(1:36)]


# COLOURS ----------------------------------------------------------------------
# This script deals with colours in sRGB values, which is a representation suitable for humans (with normal trichromatic vision).

# With all Rosaceae species ----------------------------------------------------

df_fam = subset(df, family == "Rosaceae")
#rm(df)

# converting to Lab colour space
# convert sRGB median values into Lab colour space.
Lab <- convertColor(df_fam[, c("sR_median","sG_median","sB_median")] / 255, from = "sRGB", to = "Lab") 
df_fam = cbind(df_fam, Lab)

fam_med_col = data.frame(L = with(df_fam, tapply(L, scientificName, median)),
                         a = with(df_fam, tapply(a, scientificName, median)),
                         b = with(df_fam, tapply(b, scientificName, median))
)

# converting back to sRGB colour space for display
sRGB <- convertColor(fam_med_col[,c(1:3)], from = "Lab", to = "sRGB")
fam_med_col = cbind(fam_med_col, sRGB)
names(fam_med_col)[c(4:6)] <- c("R", "G", "B")
fam_med_col$hex <- with(fam_med_col, rgb(R, G, B))
# fam_med_col has the median colour of each taxa within the family in Lab, sRGB and HEX formats


# visualise median colours
#pdf(file = "Rosaceae_sp_median_colour.pdf", width = 2.5, height = 6) 
ggplot(fam_med_col, aes(y = row.names(fam_med_col), fill = rgb(R, G, B))) +
  geom_bar(size=0.001, color = "grey20") +
  scale_y_discrete(limits = rev) +
  scale_fill_identity() +
  xlab("Median colour") +
  ylab("") +
  theme(legend.position="none", panel.grid = element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size=8), axis.ticks = element_blank(), axis.text.y = element_text(size=6), panel.background = element_rect(fill = "black"))
#dev.off()

# add colour min and max range -------------------------------------------------

df_fam$scientificName = as.factor(df_fam$scientificName)

# perform PCA on Lab values
df_min = as.data.frame(levels(df_fam$scientificName))
names(df_min)[1] = "scientificName"
df_min$scientificName = as.factor(df_min$scientificName)
df_min[c('R', 'G', 'B')] <- NA
df_max = df_min

# extract extreme RGB values of PC1
for (i in 1:length(levels(df_min$scientificName))) { 
  df_temp = subset(df_fam, scientificName == levels(df_fam$scientificName)[i])
  pca_temp = prcomp(df_temp[, c("L","a","b")], center=TRUE, scale=TRUE)
  df_temp$PC1 = pca_temp$x[,1]
  df_max$R[i] = paste(df_temp$sR_median[which.max(df_temp$PC1)]/ 255)
  df_max$G[i] = paste(df_temp$sG_median[which.max(df_temp$PC1)]/ 255)
  df_max$B[i] = paste(df_temp$sB_median[which.max(df_temp$PC1)]/ 255)
  df_min$R[i] = paste(df_temp$sR_median[which.min(df_temp$PC1)]/ 255)
  df_min$G[i] = paste(df_temp$sG_median[which.min(df_temp$PC1)]/ 255)
  df_min$B[i] = paste(df_temp$sB_median[which.min(df_temp$PC1)]/ 255)
}

df_max$type = 'max median'
df_min$type = 'min median'
df_min[,2:4] = sapply(df_min[,2:4], as.numeric)
df_max[,2:4] = sapply(df_max[,2:4], as.numeric)

# make sure the values are ordered in the same direction
df_max$type = ifelse(rowSums(df_max[,2:4]) < rowSums(df_min[,2:4]), paste('min median'), df_max$type)
df_min$type = ifelse(rowSums(df_max[,2:4]) < rowSums(df_min[,2:4]), paste('max median'), df_min$type)
fam_cols = fam_med_col[, c("R","G","B")]
fam_cols$scientificName = as.factor(rownames(fam_cols))
rownames(fam_cols) = NULL
fam_cols$type = 'median'
fam_cols = rbind(fam_cols, df_min, df_max)


# visualise min, median and max colours
#pdf(file = "Fam_colours.pdf", width = 4.3, height = 9) 
ggplot(fam_cols, aes(y = scientificName, fill = rgb(R, G, B))) +
  geom_bar(color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_identity() +
  facet_grid(.~type) +
   xlab("") +
  ylab("") +
  theme(legend.position="none", panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 8), panel.background = element_rect(fill = "black"), strip.background =element_blank(), strip.text = element_text(size = 8), panel.spacing = unit(0, "lines"))
#dev.off()

## DOMINANT COLOURS ------------------------------------------------------------

names(df_fam)

df_col_0 = df_fam[,c('image_name', 'scientificName', 'rgb_0_r', 'rgb_0_g', 'rgb_0_b', 'rgb_0_frac')]
df_col_1 = df_fam[,c('image_name', 'scientificName', 'rgb_1_r', 'rgb_1_g', 'rgb_1_b', 'rgb_1_frac')]
df_col_0$col_group = 0
df_col_1$col_group = 1

df_col_0$sum_RGB = NA
df_col_1$sum_RGB = NA

# make sure similar values are in the same class
df_col_0$sum_RGB = ifelse(rowSums(df_col_0[,3:5]) < rowSums(df_col_1[,3:5]), 0, 1)
df_col_1$sum_RGB = ifelse(rowSums(df_col_0[,3:5]) < rowSums(df_col_1[,3:5]), 1, 0)

#change column names to merge datasets
colnames(df_col_1)[3] <- paste("rgb_0_r")
colnames(df_col_1)[4] <- paste("rgb_0_g")
colnames(df_col_1)[5] <- paste("rgb_0_b")
colnames(df_col_1)[6] <- paste("rgb_0_frac")

df_col = rbind.data.frame(df_col_0, df_col_1)

colnames(df_col) = c("image_name", "scientificName", "rgb_r", "rgb_g", "rgb_b", "rgb_frac", "col_group", "sum_RGB")
str(df_col)

# convert sRGB to Lab
Lab <- convertColor(df_col[,c("rgb_r", "rgb_g", "rgb_b")]/255, from="sRGB", to="Lab")
df_col = cbind(df_col, Lab)

dom_col_sum = data.frame(
  L = with(df_col, tapply(L, list(scientificName, sum_RGB), median)),
  a = with(df_col, tapply(a, list(scientificName, sum_RGB), median)),
  b = with(df_col, tapply(b, list(scientificName, sum_RGB), median)),
  frac = with(df_col, tapply(rgb_frac, list(scientificName, sum_RGB), median))
)
dom_col_sum$scientificName = rownames(dom_col_sum)

df_dom_cols = reshape(dom_col_sum, 
                      idvar = 'scientificName',
                      varying = c(1:8) ,
                      direction = 'long')

sRGB <- convertColor(df_dom_cols[,c("L", "a", "b")], from="Lab", to="sRGB")
colnames(sRGB) = c("R", "G", "B")
df_dom_cols = cbind(df_dom_cols, sRGB)

colnames(df_dom_cols)[2] = "type"
df_dom_cols$type <- ifelse(df_dom_cols$type =="0", "colour 2", "colour 1")

# visualise dominant colours
#pdf(file = "fam_dominant_colours.pdf", width = 3.5, height = 9) 
ggplot(df_dom_cols, aes(y = scientificName, fill = rgb(R, G, B))) +
  geom_bar(color = "black") +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +
  facet_grid(.~type) +
  xlab("") +
  ylab("") +
  theme(legend.position="none", panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 8), panel.background = element_rect(fill = "black"), strip.background =element_blank(), strip.text = element_text(size = 8), panel.spacing = unit(0, "lines"))
#dev.off()

# visualise dominant colours with percentage area cover
#pdf(file = "fam_dominant_colours_frac.pdf", width = 3.5, height = 6) 
ggplot(df_dom_cols, aes(y = frac*100, x = scientificName, fill = rgb(R, G, B))) +
  geom_col() +
  scale_x_discrete(limits = rev) +
  scale_fill_identity() +
  xlab("") +
  ylab("Area cover of dominant colours (%)") +
  scale_y_continuous(limits = c(0, 101), breaks = c(0, 25, 50, 75, 100)) +
  coord_flip() +
  theme(legend.position="none", panel.grid = element_blank(), axis.text.y = element_text(size=6),  axis.text.x = element_text(size=6), axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"), strip.background =element_blank(), axis.title.x = element_text(size=6))
#dev.off()

# plot median, min, max, colour 1 and colour 2 ---------------------------------
# of Rosaceae species in a single figure ---------------------------------------

# Plot together: fam_cols and df_dom_cols

fam_all_cols = df_dom_cols[,c("scientificName", "type", "R", "G", "B")]
fam_all_cols = rbind(fam_all_cols, fam_cols)
# fam_all_cols$type = as.factor(fam_all_cols$type)
# levels(fam_all_cols$type)
fam_all_cols$type = factor(fam_all_cols$type, levels=c("max median", "median", "min median", "colour 1", "colour 2"))

# plot all colours
pdf(file = "colour_humans_Rosaceae.pdf", width = 5.5, height = 9) 
ggplot(fam_all_cols, aes(y = scientificName, fill = rgb(R, G, B))) +
  geom_bar(color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_identity() +
  facet_grid(.~type) +
  xlab("") +
  ylab("") +
  theme(legend.position="none", panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 8), panel.background = element_rect(fill = "black"), strip.background =element_blank(), strip.text = element_text(size = 8), panel.spacing = unit(0, "lines"))
dev.off()

