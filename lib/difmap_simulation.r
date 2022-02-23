# difmap_simulation.r
# Testing Diffusion Mapping on fake data
# Knights Lab - University of Minnesota
# February 2020
# usage : source('difmap_simualtion.r')

##### Set Up #####
library(diffusionMap)

##### Load Data #####
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun_long.txt", sep="\t", header=T)
dat <- dat[1:100,] # remove noise
dat <- dat[,colSums(dat) > 0]

##### Run Diffusion Mapping #####
# compute pairwise distances
bray <- dist(dat, method = "euclidean")
# diffuse computes diffusion map coordinates using pair-wise distances.
dmap <- diffuse(bray)
dcoor <- data.frame(dmap$X)
colnames(dcoor)  <- paste0("Dim", 1:ncol(dcoor)); rownames(dcoor) <- rownames(dat);
dcoor$Color <- 1:nrow(dcoor)
ggplot(dcoor, aes(x=Dim1, y=Dim2, color=factor(Color), fill=factor(Color))) + 
  geom_point(alpha=0.4, size=5, pch=21) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(title = "Diffusion Mapping",
       x = "Dimension 1",
       y = "Dimension 2") +
  guides(color="none", fill="none") +
  theme_classic() + NULL



