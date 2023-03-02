# Exploring the Y shape in PCoA
# Started with a conversation over email with Abby Johnson
# Last updated: February 2023

library(scales)
source('lib/lgd_source.r')
library(plotly)
library(vegan)

# Load Turkey Cecum dataset
dat <- readRDS("data/cecum/cecum_bray_dist.rds") # Bray-Curtis distances
meta <- read.table("data/cecum/cecum_meta.csv", sep=",", row=1, header=T)
## get rid of 7 samples removed for the paper
remove <- c("WPC2FHBD15B06C","WPC2FHBD18B07C","WPC2FHBD15B05C","WPC2FD15B04C",
            "WPC2FD15B05C","WPC2FD08B08C","WPC2FD01B10C")
remove_idx <- which(meta$SeqID %in% remove)
dat <- as.dist(as.matrix(dat)[-remove_idx,-remove_idx])
meta <- meta[-remove_idx,]
## colors and shapes
cecum_col1 <- alpha(c("#fd2beb","#8e0013","#fa4d3f","#f99314","#00ae03","#57f89f","#0e2dfc","#6033f7"),0.6)
cecum_shp1 <- c(16, 17)
ord_c <- as.factor(meta$Collection)
ord_s <- as.factor(meta$System)
cecum_col = cecum_col1[ord_c]; cecum_shp = cecum_shp1[ord_s]


# Y-shape can be created with small r here
par(mfrow=c(1,2))
for (r in c(1,0.9,0.8,0.6,0.4,0.2)) {
  lgd <- lg.dist(dat, neighborhood.radius=r)
  pc <- cmdscale(lgd, k=3, eig=T)
  plot(density(lgd), lwd=3, col="black", main=paste0("Dist Density r = ",r),
       xlab="Distance", ylab="Frequency")
  plot(pc$points, cex=2, col=cecum_col, bg=cecum_col, pch=cecum_shp,
       main=paste0("PCoA r = ",r), xlab="PC 1", ylab="PC 2",
       xlim=c(min(pc$points),max(pc$points)), ylim=c(min(pc$points),max(pc$points)))
}

# Look at the y-shape in 3D
quick_df <- data.frame(PC1 = pc$points[,1],
                       PC2 = pc$points[,2],
                       PC3 = pc$points[,3],
                       colorby = ord_c,
                       shpby = ord_s)
fig3d <- plot_ly(quick_df, x=~PC1, y=~PC2, z=~PC3, color=~colorby,
                 colors=cecum_col1, symbols=cecum_shp1)
fig3d <- fig3d %>% add_markers()
fig3d
plot(pc$eig[1:10])


# Why is this happening??? How can we avoid it?
## what about nmds with this?
nmds_out <- metaMDS(lgd, k=2, try=20, trymax=50, maxit=500)
plot(nmds_out$points, cex=2, col=cecum_col, bg=cecum_col, pch=cecum_shp,
     main="NMDS r = 0.2", xlab="NMDS 1", ylab="NMDS 2",
     xlim=c(min(nmds_out$points),max(nmds_out$points)),
     ylim=c(min(nmds_out$points),max(nmds_out$points)))
## --- also spikey but shows all in 2 dimensions this time...

## Check out qq plots
qqnorm(dat)
qqnorm(lgd)
## --- seems to follow a log-normal distribution in the spikey output
ln_data <- np.random.





