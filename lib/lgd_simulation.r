# lgd_simulation.r
# Local Gradient Distance on Simulated PCA data
# Knights Lab - University of Minnesota
# September 2019
# usage : lgd_simulation.r

##### Set Up #####
#library(optparse)
source('/project/flatiron2/suzie/detrending/fake/lgd_source.r')
library(vegan)


##### Parse Command Line #####
# option_list <- list(
#   make_option(c("-i","--input"), type="character",
#               help="Path to input distance matrix [required]."),
#   make_option(c("-n","--neighborhood_size"), type="numeric",
#               default=NULL,
#               help="Neighborhood size in which to trust distances. Smaller is better for detrending. [default: find minimum neighborhood size for connected graph]."),
#   make_option(c("-o","--output"), type="character",
#               help="Path to output distance matrix [required].")
# )
# opts <- parse_args(OptionParser(option_list=option_list), 
#                    args=commandArgs(trailing=TRUE))
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun_long.txt", sep="\t", header=T)
d <- vegdist(dat) # distance matrix from dat
neighborhood_size <- NULL


##### Helpful Functions #####
"plot.pcoa" <- function(pc, title, xlim, ylim) {
  # Note this is tailored to fake data set!
  require('ggplot2', warn.conflicts=FALSE, quietly=TRUE)
  df <- data.frame(X=pc[,1], Y=pc[,2], color=rep(1:200,2), noise=c(rep("no", 200), rep("yes", 200)))
  ggplot(df, aes(x=X, y=Y, color=factor(color), shape=factor(noise), size=factor(noise))) +
    geom_point(alpha=0.5) +
    scale_color_viridis_d(nrow(pc)/2) +
    scale_shape_manual(values=c(1, 18)) +
    scale_size_manual(values=c(15, 10)) +
    labs(title=title,
         x="PC1",
         y="PC2",
         shape="noise") +
    xlim(xlim) + ylim(ylim) +
    guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
    theme_classic() + NULL
}



##### Run PCA #####
pca <- cmdscale(vegdist(dat))
plot(pca, main = "Bray-Curtis Dissimilarity", xlab = "PC1", ylab= "PC2")
# According to PCA, the end points of this curve should be more cloesly related to one another,
#   but in reality we know these samples are at opposite ends of the OTU gradient. Can we derive
#   this gradient ourselves?


##### Local Gradient Distance #####
# Given an ordinal space, can we derive a gradient via neighbors?
cat('Calculating local gradient distance...\n')
lgd <- lg.dist(d,neighborhood.size=neighborhood_size,weighted=TRUE)
if(is.null(lgd)){
  stop("Error: graph was not connected with given neighborhood size. Try a larger neighborhood size or automatic selection.")
}
cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd <- cmdscale(lgd)
cat('Plotting PCoA of original distances...\n')
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances")
plot.pcoa(pc.d, "Original Distances", xlim=range(pc.d), ylim=range(pc.d))
cat('Plotting PCoA of transformed distances...\n')
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances")
plot.pcoa(pc.lgd, "Transformed Distances", xlim=range(pc.lgd), ylim=range(pc.lgd))



