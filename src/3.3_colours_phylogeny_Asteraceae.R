#Fede version
#4 June 2023

# devtools::install_github("jinyizju/V.PhyloMaker2")
# devtools::install_github("YuLab-SMU/ggtree")
# devtools::install_github("YuLab-SMU/treeio")
# devtools::install_github("YuLab-SMU/ggtreeExtra")
# devtools::install_github("YuLab-SMU/TDbook")
suppressPackageStartupMessages(library("V.PhyloMaker2")) # load the package
library(ggtree)
library(treeio)
library(phylobase)
library(ggtreeExtra)
library(ggnewscale)
library(factoextra)
library(phytools)
library(ape)

#setwd("write here the path to your working directory") # setting the working directory folder

df = read.csv("DiasMorph_quantitative_traits.csv")  # available at https://doi.org/10.6084/m9.figshare.21206507.v3
#names(df)

# reducing size of df. Keepings liner sRGB for analysis and sRGB for plotting.
df = df[,c("scientificName", "genus", "family", "R_median", 
           "G_median", "B_median", "R_mean", "G_mean", "B_mean", 
           "sR_median", "sG_median", "sB_median")] 

df_fam = subset(df, family == "Asteraceae")
#rm(df)

#  Prepare dataframe for V.PhyloMaker2 -----------------------------------------
sp.list = unique(df_fam[,c("scientificName", "genus", "family")])
sp.list = droplevels(sp.list)
names(sp.list)[1] = "species"

# build tree -------------------------------------------------------------------
result <- V.PhyloMaker2::phylo.maker(sp.list = sp.list, 
                                     nodes=nodes.info.1.TPL, 
                                     scenarios="S3") 
# checking tree - binary
tree = result$scenario.3
is.binary(tree)
is.binary(multi2di(tree))
tree1 <- multi2di(tree)
tree2 <- collapse.singles(tree1,root.edge=TRUE)
is.binary(tree2) # TRUE

# checking tree - ultrametric
is.ultrametric(tree2)
phy <- phangorn::nnls.tree(cophenetic(tree2), tree2, rooted = TRUE)
is.ultrametric(phy) #TRUE

# preparing labels and colours for plots ---------------------------------------

# abbreviating species name
d = sp.list[,c(1,2)]
d$temp1 = sapply(strsplit(d$species," "), `[`, 1)
d$temp1 = paste(substr(d$temp1,1,1), ".", sep = "")
d$temp2 = sapply(strsplit(d$species," "), `[`, 2)
d$temp3 = sapply(strsplit(d$species," "), `[`, 3)
d$temp4 = sapply(strsplit(d$species," "), `[`, 4)

d$label3 = paste(d$temp1, d$temp2, d$temp3, d$temp4, sep = " ")
d$label3 = sapply(strsplit(d$label3," NA"), `[`, 1)
d$species = gsub(" ", "_", d$species)

d = d[, -which(names(d) %in% c("temp1", "temp2", "temp3", "temp4"))]

colnames(d)[1] = "label"
colnames(d)[2] = "label2"
colnames(d)[3] = "label3"

rename_taxa(phy, d, label, label3) %>% write.tree
phy2 = full_join(phy, d, by = "label")
phy2

# median colours for phylo analysis and  visualisation 
# use sRGB values for tree plot
Lab <- convertColor(df_fam[,c("sR_median", "sG_median", "sB_median")] / 255, 
                    from = "sRGB", to = "Lab")
df_fam = cbind(df_fam, Lab)
sum_med_col = data.frame(R = with(df_fam, tapply(R_median, scientificName, median)),
                         G = with(df_fam, tapply(G_median, scientificName, median)),
                         B = with(df_fam, tapply(B_median, scientificName, median)),
                         L = with(df_fam, tapply(L, scientificName, median)),
                         a = with(df_fam, tapply(a, scientificName, median)),
                         b = with(df_fam, tapply(b, scientificName, median)))
sum_med_col$species = rownames(sum_med_col)

# converting back to sRGB colour space after taking the median 
sRGB <- convertColor(sum_med_col[,c("L", "a", "b")], from = "Lab", to = "sRGB")
sum_med_col = cbind(sum_med_col, sRGB)
names(sum_med_col)[c(8:10)] <- c("sR", "sG", "sB")
head(sum_med_col) #ready for plotting

# Prepare traits for phylogenetic analysis -------------------------------------

# Traits in phylogenetic order
phylo_order = as.data.frame(phy$tip.label)
phylo_order$phylo_order = as.numeric(rownames(phylo_order))
names(phylo_order)[1] = c("species")
sum_med_col$species = gsub(" ", "_", sum_med_col$species)
phylo_order = merge(phylo_order, sum_med_col, by = "species")
phylo_order = phylo_order[order(phylo_order$phylo_order),]
row.names(phylo_order) = NULL

# PCA for dimentionality reduction
pca = prcomp(phylo_order[, c("R","G","B")], center = TRUE, scale. = TRUE)
phylo_order$PC1 = pca$x[,1]
phylo_order$PC2 = pca$x[,2]

eig.val_dom <- get_eigenvalue(pca) # Eigenvalues for plot

pca_plot = ggplot(phylo_order, aes(x=PC1, y=PC2)) +
   geom_point(size = 4, shape = 16, aes(col = rgb(sR, sG, sB)), 
               position = position_nudge(x = 0.1)) +
    scale_color_identity() +
    xlab(paste("PC1  (",round(eig.val_dom[1,2], digits = 1),"%)", sep="")) +
    ylab(paste("PC2  (",round(eig.val_dom[2,2], digits = 1),"%)", sep="")) +
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 12), axis.text = element_text(size = 11))

pdf(file = "colour_phylo_pca.pdf", width = 5.5, height = 3.5)
pca_plot
dev.off()

# matching labels and PC1 scores
traits = setNames(phylo_order$PC1,
                  phylo_order$species)

# Phylogenetic signal ----------------------------------------------------------

lambda <- phylosig(phy, traits, method = "lambda", test=TRUE)
lambda

## likelihood ratio test to assess whether the phylogenetic structure fits the Brownian model. Data not shown in paper.
fitBrownian<-brownie.lite(paintSubTree(phy,phy$edge[1,1],state="1"),traits)

LR<-2*(lambda$logL-fitBrownian$logL1)
LR

P.lr<-pchisq(LR,df=1,lower.tail=FALSE)
P.lr # p<0.0001, reject the hypothesis of fitting Brownian model

# plotting tree and colours ----------------------------------------------------
# using sRGB values for colour visualisation

tr2 <- phylo4d(tree, phylo_order)

par(mar = c(2, 2, 2, 2))
p = ggtree(tr2, layout = "circular", 
           ladderize = TRUE, size=0.5) + 
    theme(legend.position="none")
#p

p2 = p %<+% d + geom_tiplab(aes(label=label3), color = "black", fontface = 3,
                             size = 1.8, h=-0.1, offset = 10, align = TRUE) 
                            #may need to adjust offset
p2

p3 = 
    p2 +
    geom_fruit(
    geom = geom_tile,
    mapping = aes(y = species, fill = rgb(sR, sG, sB)),
    #size = 5
    pwidth = 0.1,
    width = 10,
    offset = 0.12 # -0.88 (try this value for mac) offset may need to be adjusted  
  ) +
  scale_fill_identity() +
  new_scale_fill()

#p3

p4 = p3 + theme(plot.margin = margin(1.5, 1, 3, 1, "cm"))

pdf(file = "colour_phylo_tree.pdf", width = 8, height = 8)
p4
dev.off()

# plot tree with PC1 -----------------------------------------------------------

dd = phylo_order
row.names(dd) = dd$species

pr = ggtree(tree, layout = "rectangular", 
            ladderize = TRUE, size=0.6) + 
  theme(legend.position="none") #+ theme(plot.margin = margin(.3, .3, .3, .3, "cm"))

pr2 <- pr + 
  geom_facet(panel = "PC1", 
               data = dd, geom = geom_bar, fill = "grey25",
               stat = "identity", orientation = 'y',
               mapping = aes(x = PC1)) +
  theme_tree2() +
  theme(axis.text.x = element_text(size = 18), strip.text = element_blank())


pr2

# Create a data frame for the labels
d <- data.frame(.panel = c('Tree', 'PC1'), 
                lab = c("", "PC1 score"), 
                x=c(2.5,0), y=-27)

# Add text labels to the plot
pr3 <- pr2 + scale_y_continuous(limits=c(0, 195), 
                                expand=c(0,0), 
                                oob=function(x, ...) x) +
  geom_text(aes(x = x, y = y, label = lab), 
            data = d, size = 5.8, hjust = 0.5, vjust = -3.1) + 
  coord_cartesian(clip='off')  + 
  theme(plot.margin=margin(0, 10, 25, 6)) 

pr3

# Save the plot to a PDF file
pdf(file = "colour_phylo_PC1.pdf", width = 9, height = 8)
pr3
dev.off()

