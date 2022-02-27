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
plot.pc1 <- function (pc_cmd) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  pc_cmd$sample <- as.numeric(gsub("sample.", "", rownames(pc_cmd)))
  p <- ggplot(pc_cmd, aes(x=sample, y=Dim1, color=factor(color))) +
    geom_point(size = 3) +
    scale_color_manual(values = unname(colors[rownames(pc_cmd)])) +
    labs(title = "", x = "", y = "") + 
    theme_classic() + 
    theme(legend.position = 'none')
  return(p)
}

##### Loading Data #####
dat_all <- read.delim(paste0(current_dir, "data/sim_gradient/fake_rel_abun_long.txt"), sep="\t", header=T)
dat <- dat_all[1:100,] # remove noise: 1:200, all to include noise.
colors <- viridis(nrow(dat))
names(colors) <- rownames(dat)

##### Run CA on Chi-Square distances #####
# CCA function uses Chi-Sq distance and CA via SVD
ca_out <- cca(dat)
# Plot results
plot.pc1(ca_out$CA$u)