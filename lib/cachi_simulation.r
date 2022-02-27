# cachi_simulation.r
# Correspondence Analysis w/ Chi-Squared Dists on Simulation Data
# Knights Lab - University of Minnesota
# February 2022
# usage : cachi_simulation.r

##### Set Up #####
#library(optparse)
library(vegan)
library(ggplot2)
library(viridis)
library(plot.matrix)
current_dir = "~/Desktop/microbeta/"

##### Helpful functions #####
plot_dim12 <- function (pc_cmd) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  pc_cmd$sample <- as.numeric(gsub("sample.", "", 1:nrow(pc_cmd)))
  p <- ggplot(pc_cmd, aes(x=Dim1, y=Dim2, color=factor(color))) +
    geom_point(size = 3) +
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

##### Run CA on Chi-Square distances #####
# ChiSq distances
chi_d <- vegdist(dat, method = "chisq")
bray_d <- vegdist(dat, method = "bray")
jacc_d <- vegdist(dat, method = "jaccard")
# Compare distances
plt_chi_d <- as.matrix(chi_d); plt_chi_d[lower.tri(chi_d)] <- 0;
plot(plt_chi_d, col = viridis, main="Chi-Square Distances")
plt_bray_d <- as.matrix(bray_d); plt_bray_d[lower.tri(bray_d)] <- 0;
plot(plt_bray_d, col = viridis, main="Bray-Curtis Distances")
plt_jacc_d <- as.matrix(jacc_d); plt_jacc_d[lower.tri(jacc_d)] <- 0;
plot(plt_jacc_d, col = viridis, main="Jaccard Distances")
# CCA function uses Chi-Sq distance and CA via SVD
ca_out <- cca(dat)
# Plot results
plt_ca <- ca_out$CA$u
plot_dim12(plt_ca[,1:2])
