# devtools::install_github("cmartin/ggConvexHull")
library(ggplot2)
library(gridExtra)
library(factoextra)
library(ggConvexHull)

#setwd("write here the path to your working directory") # setting the working directory folder

df = read.csv("DiasMorph_v1.csv")  # available at XXXXX
# head(df)

names(df)

# reducing size of df. If working with a different version of DiasMorph, make sure it has identifying columns and colour measurements
df = df[,c(1:36)]


# COLOURS ----------------------------------------------------------------------
# This script deals with dominant colours in sRGB values, which is a representation suitable for humans (with normal trichromatic vision).


# variation of dominant colours for one taxa -----------------------------------

df_taxa = subset(df, scientificName == "Fragaria vesca")
#df_taxa = subset(df, scientificName == "Carex mucronata") #another interesting example of multicoloured seeds
#rm(df)

df_col_0 = df_taxa[,c('image_name', 'scientificName', 'rgb_0_r', 'rgb_0_g', 'rgb_0_b', 'rgb_0_frac')]
df_col_1 = df_taxa[,c('image_name', 'scientificName', 'rgb_1_r', 'rgb_1_g', 'rgb_1_b', 'rgb_1_frac')]
df_col_0$col_group = 0
df_col_1$col_group = 1

df_col_0$sum_RGB = 0
df_col_1$sum_RGB = 0

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

# we want to keep the two colours grouped so the colour combination (and not individual colour) is considered on the PCA
df_col_w = reshape(df_col, 
                    idvar = c("image_name", "scientificName"),     
                    timevar = "sum_RGB",
                    direction='wide')

# PCA --------------------------------------------------------------------------

head(df_col_w)

# perform PCA on Lab values
PCA_col_dom = prcomp(df_col_w[,c("L.0" ,"a.0", "b.0", "L.1" ,"a.1", "b.1")], center=TRUE, scale=TRUE)
df_col_w$PC1 = PCA_col_dom$x[,1]
df_col_w$PC2 = PCA_col_dom$x[,2]

# Eigenvalues
eig.val_dom <- get_eigenvalue(PCA_col_dom)

# compute the k-means clustering -----------------------------------------------

# Initialize total within sum of squares error: wss
wss <- 0
set.seed(100)
# k-means no PCA
# For 1 to 15 cluster centers
for (i in 1:15) {
  km.out <- kmeans(df_col_w[,c("L.0" ,"a.0", "b.0", "L.1" ,"a.1", "b.1")], centers = i, nstart = 1)
  # Save total within sum of squares to wss variable
  wss[i] <- km.out$tot.withinss
}

# Plot total within sum of squares vs. number of clusters
plot(1:15, wss, type = "b",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

set.seed(100)
K = kmeans(df_col_w[,c("L.0" ,"a.0", "b.0", "L.1" ,"a.1", "b.1")], 4) # for Fragaria vesca
#K = kmeans(df_col_w[,c("L.0" ,"a.0", "b.0", "L.1" ,"a.1", "b.1")], 3) # for Carex mucronata

df_col_w$kcluster = K$cluster

col_cluster = data.frame(
  PC1 = with(df_col_w, tapply(PC1, kcluster, mean)),
  PC2 = with(df_col_w, tapply(PC2, kcluster, mean)),
  kcluster = with(df_col_w, tapply(kcluster, kcluster, mean))
)

# visualise colour clusters  ---------------------------------------------------
g1 =
  ggplot(df_col_w, aes(x=PC1, y=PC2, group = kcluster)) +   
    geom_convexhull(aes(color = cluster), alpha = 0.4, fill = "grey", color = "grey38", lwd = 0.1) +
    geom_point(size=5, shape = 15, aes(col = rgb(rgb_r.0, rgb_g.0, rgb_b.0, maxColorValue = 255)), 
               position = position_nudge(x = 0.1)) +
    geom_point(size=5, shape = 15, aes(col = rgb(rgb_r.1, rgb_g.1, rgb_b.1, maxColorValue = 255)), 
               position = position_nudge(x = -0.1)) +
    geom_point(data = col_cluster, aes(x = PC1, y = PC2), 
               col = "grey28", fill = "white", shape = 21, alpha = 0.8, size = 10) +
    geom_text(data = col_cluster, aes(x = PC1, y = PC2, label = kcluster), 
              size = 6, fontface = "bold") +
    scale_color_identity() +
    xlab(paste("PC1  (",round(eig.val_dom[1,2], digits = 1),"%)", sep="")) +
    ylab(paste("PC2  (",round(eig.val_dom[2,2], digits = 1),"%)", sep="")) +
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 12), axis.text = element_text(size = 11))

g1

# calculate median dominant colours for each cluster ---------------------------
df_col_w$cluster = K$cluster
k_col_dom_wide = as.data.frame(K$centers)

col_freq = data.frame(
            cluster = with(df_col_w, tapply(cluster, cluster, mean)),
                      number_samples = with(df_col_w, tapply(cluster, cluster, length)),
                      freq_dom.0 = with(df_col_w, tapply(rgb_frac.0, cluster, median)),
                      freq_dom.1 = with(df_col_w, tapply(rgb_frac.1, cluster, median)))

k_col_dom_wide = cbind(k_col_dom_wide, col_freq)

k_col_dom = reshape(k_col_dom_wide,
                    varying = c("L.0", "a.0", "b.0", "L.1", "a.1", "b.1",
                                "freq_dom.0", "freq_dom.1"),
                    direction = 'long')


sRGB <- convertColor(k_col_dom[,c("L", "a", "b")], from ="Lab", to="sRGB")
k_col_dom = cbind(k_col_dom, sRGB)  
colnames(k_col_dom)[which(names(k_col_dom) == "time")] <- "col_group"
colnames(k_col_dom)[which(names(k_col_dom) == "1")] <- "R"
colnames(k_col_dom)[which(names(k_col_dom) == "2")] <- "G"
colnames(k_col_dom)[which(names(k_col_dom) == "3")] <- "B"
k_col_dom$hex <- with(k_col_dom, rgb(R, G, B))

k_col_dom$frac_samples = round(with(k_col_dom, 
                                 tapply(number_samples, cluster, sum)) / 
                                 sum(k_col_dom$number_samples) * 100, 1)

# visualise median colour for each cluster -------------------------------------
g2 =
  ggplot(k_col_dom, aes(y = freq_dom * 100, x = cluster, fill = rgb(R, G, B))) +
    geom_col() +
    geom_text(aes(x = cluster, y = 100, label = paste(format(frac_samples, nsmall = 1),"%")), 
              size = 3.5, hjust = -0.1) +
    scale_fill_identity() +
    xlab("Cluster identity") +
    ylab("Area cover of\ndominant colours (%)") +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 25, 50, 75, 100)) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank(), 
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), 
          axis.text = element_text(size = 11)) +
    coord_flip()

g2

# save PDF ---------------------------------------------------------------------
pdf(file = "colour_dominant_Fragaria_vesca.pdf", width = 9, height = 4) 
grid.arrange(arrangeGrob(g1, g2, as.table= TRUE, ncol=2, nrow=1, heights=c(5), widths=c(6.5, 2.8)))
dev.off()
