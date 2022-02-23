# isomap_simulation.r
# IsoMap on Simulation Data
# Knights Lab - University of Minnesota
# July 2019
# usage : isomap_simulation.r -i input_table.txt

##### Set Up #####
#library(optparse)
library(vegan)
library(ggplot2)


##### Parse Command Line #####
# option_list <- list(make_option(c("-i", "--input_table"), type="character",
#                                 help="Path to input file. Expects otu table in tab-delimited txt file format.") )
# opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun.txt", sep="\t", header=T)


##### Helpful Functions #####
colorfunc <- colorRampPalette(c("white", "navy"))
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}


##### Run Isomap #####
distances <- dist(dat, method = "euclidean")
iso_out <- isomap(distances, ndim=2, epsilon=0.3)
ggplot(data.frame(X=iso_out$points[,1],Y=iso_out$points[,2], color=rep(1:200,2), noise=c(rep("no",200), rep("yes",200))),
       aes(x=X, y=Y, color=factor(color), shape=factor(noise), size=factor(noise))) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(nrow(iso_out$points)/2) +
  scale_shape_manual(values=c(1, 18)) +
  scale_size_manual(values=c(15, 10)) +
  labs(title="IsoMap",
       x=paste0("Dimension 1 [", calc.perc.var(iso_out$eig,1), "%]"),
       y=paste0("Dimension 2 [", calc.perc.var(iso_out$eig,2), "%]"),
       shape="noise") +
  guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
  theme_classic() + NULL

