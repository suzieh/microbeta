# cachi_simulation.r
# Correspondence Analysis w/ Chi-Squared Dists on Simulation Data
# Knights Lab - University of Minnesota
# February 2022
# usage : cachi_simulation.r

##### Set Up #####
#library(optparse)
library(vegan)
library(ggplot2)
library(ggpubr)
current_dir = "~/Desktop/microbeta"

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
dat_all <- read.delim("fake_rel_abun.txt", sep="\t", header=T)
dat <- dat_all[1:200,] # remove noise for now
colors <- viridis(nrow(dat))
names(colors) <- rownames(dat)


