library(Momocs)
library(ggplot2)
library(sf)

#setwd("XXX") # insert working directory

df = read.csv("DiasMorph_v1.csv")  # available at XXXXX
# head(df)

# selecting only Carex spp.
df = subset(df, genus == "Carex")

# randomly selecting only taxa with more than 30 seeds
number.samples <- as.data.frame(tapply(df$scientificName, df$scientificName, length))
colnames(number.samples)[1] <- 'count'
number.samples$scientificName <- row.names(number.samples)
not_enough = number.samples[number.samples$count < 30, ]
Slevels = levels(as.factor(not_enough$scientificName))
df = df[! (df$scientificName %in% Slevels), ]

# Number of taxa to be included in the analysis
length(levels(as.factor(df$scientificName)))

# Randomly downsampling taxa with more than 30 seeds to obtain a balanced dataset
FUN <- function(x, n) {
  if (length(x) <= n) return(x)
  x[x %in% sample(x, n)]
}

set.seed(8)
df <- df[unlist(lapply(split(1:nrow(df), df$scientificName),
                       FUN, n = 30)), ]
# final dataset


# format dataset for Momocs ----------------------------------------------------
names(df)

# converting df into a list of matrices of (x,y) coordinates
coord_list = list()
for (i in 1:length(df[,1])){
  coord_list[[i]] = t(array(
    as.matrix(df[i,c(37:136)]), # select columns from x_0 to y_49
    dim = c(2,50)
  ))
  names(coord_list)[i] <- paste(df[i,"image_name"])
}

# data.frame of variables specifying the grouping structure
fac=list()
fac$sp <- df$scientificName
fac$sample <- df$image_name
fac$sp_code <- df$group

# import coordinates into momocs 
out_df = Out(coord_list, fac=fac)

# visualise all shapes by species
#panel(out_align1, cols=as.factor(fac$sp))


# aligning shapes --------------------------------------------------------------

#  perform a 90-degree rotation in shape outline
out_align <- out_df %>%
  coo_rotate(pi/2) 

## performing Elliptical Fourier Transforms followed by PCA on the coefficients
pca_align <- out_align %>% efourier(norm=FALSE) %>% PCA(scale. = FALSE, center = TRUE)
plot_PCA(pca_align, axes = c(1, 2), labelpoints = TRUE)  #ploting labels to identify any outliers

# removing outlier to improve outcome of the alignment
out_align <- filter(out_align, !sample %in% c("img_0208_6","img_0208_40"))

# use PCA to separate seeds with mirrored alignment
pca_align <- out_align %>% efourier(norm=FALSE) %>% PCA(scale. = FALSE, center = TRUE)
plot_PCA(pca_align, axes = c(1, 2)) 

# PC2 is separating seeds with apex on the left from those with apex on the right.

#Split dataset into outlines with PC2 >= 0 from those with PC2 < 0
out_alignA <- filter(out_align,pca_align$x[,2]>= 0)
#stack(out_alignA) # visualise alignment of outlines with PC2 >= 0

out_alignB <- filter(out_align,pca_align$x[,2]< 0)
#stack(out_alignB) # now compare with alignment of outlines with PC2 < 0

# Apply a 180-degree rotation to all oulines with PC2<0 to ensure that all seeds have the apex facing the same direction.
out_alignB <- out_alignB %>% coo_flipy() %>% coo_rev() # mirrors shape horizontally

# combine dataset
out_fin <- combine(out_alignA,out_alignB)
#panel(out_fin, cols=as.factor(out_fin$fac$sp)) #visualise outlines after alignment

# performing Elliptical Fourier Transforms followed by PCA on the coefficients
res_pca = out_fin %>%
  efourier(norm=FALSE) %>%
  PCA(scale. = FALSE, center = TRUE)

par(mfrow = c(1,1))
plot_PCA(res_pca, axes = c(1, 2), chull = FALSE)


# hierarchical clustering ------------------------------------------------------

k = 8 #specify number of clusters 

par(mfrow=c(1,1))
par(mar=c(3,3,3,3))
fit_hc = hclust(dist(res_pca$x[,c(1:3)]), method = "ward.D2")
plot(fit_hc, labels = FALSE)
rect.hclust(fit_hc, k = k)

# setting k clusters as output
hc_cl = factor(cutree(fit_hc, k = k))  
plot(res_pca$x[,c(1:2)], col = as.vector(hc_cl)) #visualise the clusters

out_fin[[2]]$hc <- as.factor(hc_cl)
fac$hc <- as.factor(hc_cl)


# find the representative shape of the cluster----------------------------------

# find cluster center
hcluster_center = aggregate(res_pca$x[,c(1:3)], list(cluster = hc_cl), mean)

# find nearest point to center
pca_x=as.data.frame(res_pca$x)
pca_x$name = row.names(pca_x)


rep_shape = cbind(
  hcluster_center,
  nearest_point = pca_x[
    st_nearest_feature(
      st_as_sf(hcluster_center, coords = c("PC1", "PC2")), 
      st_as_sf(pca_x, coords = c("PC1", "PC2"))
    ),
  ]$name
)


# plot the shapes closest to the center of each cluster ------------------------

col_val = c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#58C9FB", "#B4485C", "#75326B", "#37778C")

# rotate shapes so they are plotted in upright position
plot_shapes <- out_fin %>%
  coo_rotate(3*pi/2) 

# get outline of representative shapes from plot_shapes
rep_list = list()
for (i in 1:length(rep_shape$nearest_point)){
  rep_list[[i]] = filter(plot_shapes, sample == rep_shape[i,5])
}


pdf(file = "cluster_representatives.pdf", width = 1.5, height = 12) 

par(pty="s")
par(mfrow=c(1,1))
len = length(rep_shape$cluster)
layout(matrix(c(1: len), len, 1, byrow = TRUE)
)

for (i in 1:length(rep_shape$nearest_point)){
  coo_plot(coo_smooth(rep_list[[i]]$coo, 1), points = FALSE, xy.axis=FALSE 
           , col= col_val[i]
           , ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5)
  )
}

dev.off()

# panels to visualise all shapes within each cluster ---------------------------

dev.off()
pdf(file = "shape_cluster_representatives.pdf", width = 15, height = 8) 

par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1))
for (i in 1:length(rep_shape$cluster)){
  p = plot_shapes  %>%
    filter(plot_shapes[[2]]$hc == i)
  panel(p, #cols = as.factor(p[[2]]$sp)
        c(30,50),
        col = col_val[i]
  )
}

dev.off()


# a nicer PCA plot in ggplot2 --------------------------------------------------

to_plot = as.data.frame(res_pca$x[,c(1:2)])
cols = as.factor(hc_cl)

pdf(file = "shape_pca_clusters.pdf", width = 6.5, height = 4) 

ggplot() +
  geom_point(data = to_plot, aes(x = PC1, y = PC2, color = cols), 
             size = 2, alpha = 0.3) +
  geom_point(data = rep_shape, aes(x = PC1, y = PC2, color = levels(cols)),
             size = 7, shape = 4, stroke = 2.5) +
  scale_color_manual(name=" ", values=c(col_val), labels = levels(cols)) +
  ylab(paste("PC2  (",round(res_pca$eig[2]*100, digits = 1),"%)", sep="")) +
  xlab(paste("PC1  (",round(res_pca$eig[1]*100, digits = 1),"%)", sep="")) +
  theme_bw(base_size=12) +
  theme(legend.position="none")

dev.off()

# shape breadth ---------------------------------------------------------------

levels(as.factor(res_pca$fac$sp_code))
examples = filter(res_pca, res_pca$fac$sp == "Carex pauciflora" |
                    res_pca$fac$sp == "Carex pulicaris" |
                    res_pca$fac$sp == "Carex atrata subsp. aterrima" |
                    res_pca$fac$sp == "Carex filiformis" |
                    res_pca$fac$sp == "Carex divulsa" |
                    res_pca$fac$sp == "Carex viridula subsp. oedocarpa"
) 

levels(as.factor(examples$fac$sp))
sp_example = as.factor(examples$fac$sp)

ex_df = as.data.frame(examples$x[,1:2])
ex_df$sp = examples$fac$sp
ex_df$sp_code = examples$fac$sp_code

sum_ex_df = data.frame(PC1 = with(ex_df, tapply(PC1, sp, mean)),
                       PC2 = with(ex_df, tapply(PC2, sp, mean))
)

sum_ex_df$sp = rownames(sum_ex_df)

pdf(file = "shape_breadth.pdf", width = 6.5, height = 4) 

ggplot() +
  geom_point(data = to_plot, aes(x = PC1, y = PC2, color = cols), 
             size = 2, alpha = 0.5) +
  scale_color_manual(name=" ", values=col_val, labels = levels(cols)) +
  stat_ellipse(data = as.data.frame(examples$x), 
               aes(x = PC1, y = PC2, group = sp_example), 
               size = 0.4, col = "grey48", level = 0.95) +
  geom_text(data = sum_ex_df[-1,], aes(x = PC1, y = PC2, label = sp), 
            size = 3.5) +
  geom_text(data = sum_ex_df[1,], aes(x = PC1, y = PC2, label = sp), 
            size = 3.5, vjust = -3.5, hjust = 0.6) +
  ylim(-0.15, 0.10) +
  ylab(paste("PC2  (",round(res_pca$eig[2]*100, digits = 1),"%)", sep="")) +
  xlab(paste("PC1  (",round(res_pca$eig[1]*100, digits = 1),"%)", sep="")) +
  theme_bw(base_size=12) +
  theme(legend.position="none")

dev.off()
