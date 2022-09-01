# Playing around with objective functions

# Pull in cecum data
distobj <- readRDS("data/cecum/cecum_bray_dist.rds")
cecum <- read.csv("data/cecum/cecum_meta.csv", header=T)
head(cecum)

# Compute LGD
source("lib/lgd_source.r")
lgdobj <- as.dist(lg.dist(distobj, neighborhood.radius = 0.88))

# Play around with histograms
hist(distobj, breaks=30, main="Histogram of Distances", xlab="Distances")
hist(lgdobj, breaks=20, main="Histogram of LGD-adjusted", xlab="LGD Distances")
## basically, we see the density at 1.0 thin out after LGD is applied.
library(ggplot2)
ggplot(data.frame(d=c(distobj)), aes(d)) +
  geom_density(adjust=1)
ggplot(data.frame(d=c(lgdobj)), aes(d)) +
  geom_density(adjust=1)


# Testing Idea: Peak at end must fall below earlier peak(s)?
# Find peaks (know the order and proximity to maximum distance)
library(pracma)

## Visualize with Original Distances
den_do <- density(distobj)
p_do <- findpeaks(den_do$y)
plot(den_do, main="Density of Original Dists", xlab="Original Distances")
abline(v=den_do$x[p_do[,2]], col="blue")

## Visualize with Adjusted Distances
den_lgd <- density(lgdobj)
p_lgd <- findpeaks(den_lgd$y)
plot(den_lgd, main="Density of LGD Dists", xlab="LGD-adjusted distances")
abline(v=den_lgd$x[p_lgd[,2]], col="purple")

## Objective function (to minimize): Distance between peaks
##   Should measure distance between the last peak, and the maximum other peak found in density
## What if the distribution is not multi-modal???
library(stringr)
library(viridis)
obj_peak <- function (dist_obj, r) {
  den <- density(dist_obj)
  peaks <- findpeaks(den$y)
  npeaks <- nrow(peaks)
  comp_peak <- which.max(peaks[-npeaks,1])
  dif_peaks <- peaks[npeaks,1] - peaks[comp_peak,1] # last - other max
  # plotting
  png(paste0("results/peaks_gif/r_",round(r,4),".png"), width=900, height=500)
  par(mfrow=c(1,2))
  par(mar=c(5, 4, 4, 2), xpd=F)
  plot(den, main=paste0("Density of Distances: r=", round(r,4)), xlab="distances", lwd=2)
  abline(v=den$x[peaks[,2]], col="blue", lwd=2)
  mtext(paste0("diff: ",round(dif_peaks,3)), col="blue", side=3, adj=1)
  par(mar=c(5, 4, 4, 6), xpd=T)
  plot(cmdscale(dist_obj,k=2), col=alpha(viridis(8)[as.factor(cecum$Collection)],0.7), cex=2,
       pch=c(17,16)[as.factor(cecum$System)], main="PCoA", xlab="PC1", ylab="PC2")
  legend("topright", inset=c(-0.2,0), legend = c("Pen", "Hatch", levels(as.factor(cecum$Collection))),
         pch=c(16,17,rep(16,8)), col=c("black", "black", viridis(8)))
  dev.off()
  # returning error value
  if (length(dif_peaks) == 0) {
    return(Inf)
  } else if (dif_peaks > 0) {
    return(dif_peaks)
  } else {
    return(abs(dif_peaks))
  }
}

## Attempt to find best lgd with trying 100 radii b/w minimum and maximum distances
minr <- min(distobj)
maxr <- max(distobj)
bestr <- maxr
bestdist <- distobj
besterr <- Inf

# !!!!! testing source functions are working !!!!!:
# my_seq <- seq(maxr, minr, length.out=25)
# source("lib/lgd_source.r")
# my_dist <- as.dist(lg.dist(distobj, neighborhood.radius = my_seq[12]))

## Run on multiple r values to see what happens
for (r in seq(maxr, minr, length.out=25)) {
  print(r)
  tmp_dist <- as.dist(lg.dist(distobj, neighborhood.radius = r))
  # err <- obj_peak(tmp_dist, r)
  # if (err < besterr) {
  #   bestr <- r
  #   bestdist <- tmp_dist
  #   besterr <- err
  # }
}
bestr
besterr

## Make a gif from images & dif peaks calculated
library(magick)
## list file names (in order) and read in
imgs <- sapply(round(my_seq, 4), function (r) {paste0("results/peaks_gif/r_",r,".png")})
img_list <- lapply(imgs, image_read)
## join the images together
img_joined <- image_join(img_list)
## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 1)
## view animated image
img_animated
## save to disk
image_write(image = img_animated, path = "results/peaks_gif/my_peaks_gif_cecum.gif")

## Make a gif of all the igraphs too
imgs <- sapply(round(my_seq, 4), function (r) {paste0("results/igraphs_gif/r_",r,".png")})
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
img_animated
image_write(image = img_animated, path = "results/igraphs_gif/my_igraphs_gif_cecum.gif")



# Another metric: Comparing Community & PC distances
## Plotting community distances vs. pc distances (GIF again)
library(vegan)
blue_pal <- colorRampPalette(c('deepskyblue','darkblue'))
plot_comm_pc <- function (comm, pc, r) {
  pc_dist <- vegdist(pc, method = "euclidean")
  my_colors <- blue_pal(20)[as.numeric(cut(comm, breaks = 20))] # color on a gradient by community distances
  my_cor <- cor.test(comm, pc_dist, method="pearson")
  adjx <- max(comm)/10
  adjy <- max(pc_dist)/10
  plot(comm, pc_dist, col=alpha(my_colors,0.7), cex=2, pch=16,
            main=paste0("Community vs PCoA Distances: r = ", round(r, 2)),
            xlab="Community distance", ylab="PCoA distance")
  text(x = min(comm)+adjx, y = max(pc_dist)-adjy, labels = paste0("cor: ", round(my_cor$estimate, 2)), xpd = NA)
  text(x = min(comm)+adjx, y = max(pc_dist)-(2*adjy), labels = paste0("p-val: ", round(my_cor$p.value, 2)), xpd = NA)
}

## For loop to create plots for each r value
for (r in seq(maxr, minr, length.out=25)) {
  my_lgd <- as.dist(lg.dist(distobj, neighborhood.radius = r))
  my_pc <- cmdscale(my_lgd, k=3)
  png(paste0("results/dist_compare_gif/r_",round(r,4),".png"), width=500, height=500)
  plot_comm_pc(my_lgd, my_pc, r)
  dev.off()
}

## Create gif
imgs <- sapply(round(seq(maxr, minr, length.out=25), 4), function (r) {paste0("results/dist_compare_gif/r_",r,".png")})
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
img_animated
image_write(image = img_animated, path = "results/dist_compare_gif/my_distcompare_gif_cecum.gif")



# COMBINATION OBJ. FUNCTION: Maximize r and maximize correlation b/w distances
## Find the cross between r and correlation (see slide for visual representation)
##    minimize the difference between these (this will the crossing point)
## --- Set Up
minr <- min(distobj)
maxr <- max(distobj)
bestr <- maxr
bestdist <- distobj
besterr <- Inf
## --- Loop
for (r in seq(maxr, minr, length.out=25)) {
  my_lgd <- as.dist(lg.dist(distobj, neighborhood.radius = r))
  my_pc <- vegdist(cmdscale(my_lgd, k=3), method="euclidean")
  my_corr <- cor.test(my_lgd, my_pc, method="pearson")
  opt_out <- abs(r - my_corr$estimate[1])
  if (opt_out < besterr) {
    bestdist <- my_lgd
    besterr <- opt_out
    bestr <- r
  }
}
besterr
bestr
##     note that this took a while... can we speed it up??? Maybe the value compute within the LGD function.
## --- Visualize the Before & After:
##     4 plots for the original distance object
# general set up for plots
r = 1
pc <- cmdscale(distobj, k=3)
pc_dist <- vegdist(pc, method="euclidean")
adjx <- max(distobj)/10
adjy <- max(pc_dist)/10
my_cor <- cor.test(distobj, pc_dist, method="pearson")
my_colors <- blue_pal(20)[as.numeric(cut(distobj, breaks = 20))]
d_graph <- lg.graph(distobj, neighborhood.size=1, use.r=T, weighted=T)
den <- density(distobj)
par(mfrow=c(2,2))
# comm dist plot
plot(distobj, pc_dist, col=alpha(my_colors,0.7), cex=2, pch=16,
     main=paste0("Community vs PCoA Distances: r = ", round(r, 2)),
     xlab="Community distance", ylab="PCoA distance")
text(x = min(distobj)+adjx, y = max(pc_dist)-adjy, labels = paste0("cor: ", round(my_cor$estimate, 2)), xpd = NA)
text(x = min(distobj)+adjx, y = max(pc_dist)-(2*adjy), labels = paste0("p-val: ", round(my_cor$p.value, 2)), xpd = NA)
# igraph plot
plot(d_graph, vertex.label=NA, vertex.color="blue", vertex.size=8,
     main=paste0("r = ", round(r,4)))
# dist distribution plot
plot(den, main=paste0("Density of Distances: r = ", round(r,4)), xlab="distances", lwd=3, col="blue")
# pcoa plot
par(mar=c(5, 4, 4, 6), xpd=T)
plot(pc, col=alpha(viridis(8)[as.factor(cecum$Collection)],0.7), cex=2,
     pch=c(17,16)[as.factor(cecum$System)], main="PCoA", xlab="PC1", ylab="PC2")
legend("topright", inset=c(-0.43,0), legend = c("Pen", "Hatch", levels(as.factor(cecum$Collection))),
       pch=c(16,17,rep(16,8)), col=c("black", "black", viridis(8)), cex=1.3,
       y.intersp=0.3, x.intersp=0.5, xjust=0, bty = "n")
par(mfrow=c(1,1), mar=c(5, 4, 4, 2), xpd=F)
##     4 plots for "bestdist" object (using bestr)
# general set up for plots
r = bestr
pc <- cmdscale(bestdist, k=3)
pc_dist <- vegdist(pc, method="euclidean")
adjx <- max(bestdist)/10
adjy <- max(pc_dist)/10
my_cor <- cor.test(bestdist, pc_dist, method="pearson")
my_colors <- blue_pal(20)[as.numeric(cut(bestdist, breaks = 20))]
d_graph <- lg.graph(bestdist, neighborhood.size=1, use.r=T, weighted=T)
den <- density(bestdist)
par(mfrow=c(2,2))
# comm dist plot
plot(bestdist, pc_dist, col=alpha(my_colors,0.7), cex=2, pch=16,
     main=paste0("Community vs PCoA Distances: r = ", round(r, 2)),
     xlab="Community distance", ylab="PCoA distance")
text(x = min(bestdist)+adjx, y = max(pc_dist)-adjy, labels = paste0("cor: ", round(my_cor$estimate, 2)), xpd = NA)
text(x = min(bestdist)+adjx, y = max(pc_dist)-(2*adjy), labels = paste0("p-val: ", round(my_cor$p.value, 2)), xpd = NA)
# igraph plot
plot(d_graph, vertex.label=NA, vertex.color="blue", vertex.size=8,
     main=paste0("r = ", round(r,4)))
# dist distribution plot
plot(den, main=paste0("Density of Distances: r = ", round(r,4)), xlab="distances", lwd=3, col="blue")
# pcoa plot
par(mar=c(5, 4, 4, 6), xpd=T)
plot(pc, col=alpha(viridis(8)[as.factor(cecum$Collection)],0.7), cex=2,
     pch=c(17,16)[as.factor(cecum$System)], main="PCoA", xlab="PC1", ylab="PC2")
legend("topright", inset=c(-0.43,0), legend = c("Pen", "Hatch", levels(as.factor(cecum$Collection))),
       pch=c(16,17,rep(16,8)), col=c("black", "black", viridis(8)), cex=1.3,
       y.intersp=0.3, x.intersp=0.5, xjust=0, bty = "n")
par(mfrow=c(1,1), mar=c(5, 4, 4, 2), xpd=F)



