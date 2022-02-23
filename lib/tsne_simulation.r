# tsne_simulation.r
# T-SNE on Simulation Data
# Knights Lab - University of Minnesota
# July 2019
# usage : tsne_simulation.r -i input_table.txt

##### Set Up #####
#library(optparse)
library(plotly)
library(Rtsne)


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


##### Run T-SNE #####
# 2 Dimensions
tsne_out <- Rtsne(dat, dims=2, perplexity=15, verbose=T, max_iter=1000)
tsne_points <- as.data.frame(tsne_out$Y); colnames(tsne_points) <- c("Dim1", "Dim2");
tsne_points$color <- rep(1:(nrow(tsne_points)/2),2)
tsne_points$noise <- c(rep("no",(nrow(tsne_points)/2)), rep("yes",(nrow(tsne_points)/2)))
ggplot(tsne_points, aes(x=Dim1, y=Dim2, color=factor(color), shape=factor(noise), size=factor(noise))) +
  geom_point(alpha=0.5) +
  scale_color_viridis_d(nrow(tsne_points)/2) +
  scale_shape_manual(values=c(1, 18)) +
  scale_size_manual(values=c(15, 10)) +
  labs(title="t-SNE",
       x="Dimension 1",
       y="Dimension 2",
       shape="noise") +
  guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
  theme_classic() + NULL
# What happens if we expand to 3 dimensions? (removed noise)
tsne_3 <- Rtsne(dat, dims=3, perplexity=15, verbose=T, max_iter=1000)
tsne_points3 <- as.data.frame(tsne_3$Y); colnames(tsne_points3) <- c("Dim1", "Dim2", "Dim3");
tsne_points3$color <- factor(rep(1:(nrow(tsne_points)/2),2))
tsne_points3$noise <- factor(c(rep("no",(nrow(tsne_points)/2)), rep("yes",(nrow(tsne_points)/2))))
plot_ly(tsne_points3, x=~Dim1, y=~Dim2, z=~Dim3, mode="markers",
        marker=list(color = ~color, colorscale = 'Viridis', showscale = FALSE, size = 12)) %>%
  add_markers()

