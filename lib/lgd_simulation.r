# lgd_simulation.r
# Local Gradient Distance on Simulated PCA data
# Knights Lab - University of Minnesota
# September 2019
# usage : lgd_simulation.r

##### Set Up #####
#library(optparse)
source('~/Desktop/microbeta/lib/lgd_source.r')
library(vegan)
library(ggplot2)
library(viridis)
current_dir = "~/Desktop/microbeta/"


##### Helpful functions #####
plot_dim12 <- function (pc_cmd, title, flip_dim1=F, flip_dim2=F) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  pc_cmd$sample <- as.numeric(gsub("sample.", "", 1:nrow(pc_cmd)))
  if (flip_dim1) { pc_cmd$Dim1 <- -pc_cmd$Dim1 }
  if (flip_dim2) { pc_cmd$Dim2 <- -pc_cmd$Dim2 }
  p <- ggplot(pc_cmd, aes(x=Dim1, y=Dim2, color=factor(color))) +
    geom_point(size = 6, alpha=0.6) +
    scale_color_manual(values = unname(colors[rownames(pc_cmd)])) +
    labs(title = title, x = "Dimension 1", y = "Dimension 2") + 
    theme_classic() + 
    lims(y = c(-2.4,2.4)) +
    theme(legend.position = 'none',
          plot.title = element_text(size=18, face="bold", hjust = 0.5))
  return(p)
}


##### Loading Data #####
dat_all <- read.delim(paste0(current_dir, "data/sim_gradient/fake_rel_abun_long.txt"), sep="\t", header=T)
dat <- dat_all[c(rbind(seq(1,50,1),seq(51,100,1))),] # rearranging to have noise next to matching "samples".
colors <- viridis(nrow(dat))
names(colors) <- rownames(dat)


##### Distance & Distance Comparisons (density plot) #####
d <- vegdist(dat, method="bray")
plot(density(d), main="Density of Original Distances", col="blue", lwd = 3)
lgd <- lg.dist(d,neighborhood.radius=0.6,weighted=TRUE)
plot(density(lgd), main="Density of LGD Distances", col="purple", lwd = 3)

##### LGD #####
# PCoA of original distances
pc.d <- cmdscale(d)
# PCoA of transformed distances
pc.lgd <- cmdscale(lgd)
# Plot original PCoA
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances")
plot_dim12(pc.d, "Original Distances", flip_dim1 = T)
# Plot transformed distances PCoA
plot(pc.lgd, xlim=range(pc.d), ylim=range(pc.d), main="Transformed Distances")
plot_dim12(pc.lgd, "Transformed Distances")








##################################################
################ OLD CODE BELOW ##################
##################################################

##### Helpful Functions #####
# "plot.pcoa" <- function(pc, title, xlim, ylim) {
#   # Note this is tailored to fake data set!
#   require('ggplot2', warn.conflicts=FALSE, quietly=TRUE)
#   df <- data.frame(X=pc[,1], Y=pc[,2], color=rep(1:200,2), noise=c(rep("no", 200), rep("yes", 200)))
#   ggplot(df, aes(x=X, y=Y, color=factor(color), shape=factor(noise), size=factor(noise))) +
#     geom_point(alpha=0.5) +
#     scale_color_viridis_d(nrow(pc)/2) +
#     scale_shape_manual(values=c(1, 18)) +
#     scale_size_manual(values=c(15, 10)) +
#     labs(title=title,
#          x="PC1",
#          y="PC2",
#          shape="noise") +
#     xlim(xlim) + ylim(ylim) +
#     guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
#     theme_classic() + NULL
# }
# 
# 
# 
# ##### Run PCA #####
# pca <- cmdscale(vegdist(dat))
# plot(pca, main = "Bray-Curtis Dissimilarity", xlab = "PC1", ylab= "PC2")
# # According to PCA, the end points of this curve should be more cloesly related to one another,
# #   but in reality we know these samples are at opposite ends of the OTU gradient. Can we derive
# #   this gradient ourselves?
# 
# 
# ##### Local Gradient Distance #####
# # Given an ordinal space, can we derive a gradient via neighbors?
# cat('Calculating local gradient distance...\n')
# lgd <- lg.dist(d,neighborhood.radius = ,weighted=TRUE)
# if(is.null(lgd)){
#   stop("Error: graph was not connected with given neighborhood size. Try a larger neighborhood size or automatic selection.")
# }
# cat('Calculating PCoA of original distances...\n')
# pc.d <- cmdscale(d)
# cat('Calculating PCoA of transformed distances...\n')
# pc.lgd <- cmdscale(lgd)
# cat('Plotting PCoA of original distances...\n')
# plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances")
# plot.pcoa(pc.d, "Original Distances", xlim=range(pc.d), ylim=range(pc.d))
# cat('Plotting PCoA of transformed distances...\n')
# plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances")
# plot.pcoa(pc.lgd, "Transformed Distances", xlim=range(pc.lgd), ylim=range(pc.lgd))



