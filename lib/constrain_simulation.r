# constrain_simulaation.r
# Correspondence Analysis w/ Chi-Squared Dists on Simulation Data
# Knights Lab - University of Minnesota
# February 2022
# usage : cachi_simulation.r

##### Set Up #####
#library(optparse)
library(vegan)
library(ggplot2)
library(viridis)
current_dir = "~/Desktop/microbeta/"

##### Helpful functions #####
plot_dim12 <- function (pc_cmd) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  pc_cmd$sample <- as.numeric(gsub("sample.", "", 1:nrow(pc_cmd)))
  #pc_cmd$Dim1 <- -pc_cmd$Dim1
  #pc_cmd$Dim2 <- -pc_cmd$Dim2
  p <- ggplot(pc_cmd, aes(x=Dim1, y=Dim2, color=factor(color))) +
    geom_point(size = 4) +
    scale_color_manual(values = unname(colors[rownames(pc_cmd)])) +
    labs(title = "", x = "Dimension 1", y = "Dimension 2") + 
    theme_classic() + 
    theme(legend.position = 'none')
  return(p)
}

##### Loading Data #####
dat_all <- read.delim(paste0(current_dir, "data/sim_gradient/fake_rel_abun_long.txt"), sep="\t", header=T)
dat <- dat_all[c(rbind(seq(1,50,1),seq(51,100,1))),] # rearranging to have noise next to matching "samples".
colors <- viridis(nrow(dat))
names(colors) <- rownames(dat)

##### Run RDA #####
# Create a metadata variable(s)
gradient_loc <- rep(1:50,each = 2)
noise <- rep(0:1, times = 50)
# RDA computation
rda_out <- rda(dat ~ gradient_loc, scale = T) # can also add + noise
# Biplot of RDA results
rda_pts <- scores(rda_out)$sites
colnames(rda_pts) <- c("Dim1", "Dim2")
rda_arrows <- rda_out$CCA$envcentre
p <- plot_dim12(rda_pts)
p <- p + annotate("text", x = 2, y = 0.2, label = "gradient", col = "blue") +
  annotate("segment", x = 0, xend = 2, y = 0, yend = 0, col = "blue", arrow = arrow(ends = "last",length = unit(.2,"cm")))
  #annotate("text", x = 0, y = 7, label = "noise", col = "blue") +
  #annotate("segment", x = 0, xend = 0, y = 0, yend = 5, col = "blue", arrow = arrow(ends = "last",length = unit(.2,"cm")))
p

##### Run DCA #####