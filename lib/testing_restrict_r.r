# testing_restrict_r.r
# Restricting radius in LGD algorithm
# Suzie Hoops
# Created March 2023

library(scales)
set.seed(125)

# TEMP: TO-DO
# - add a restriction to radius r based on average degree


# Scratch stuff (while working on LGD source)
# Soil dataset
soil_norm <- read.table("data/soil/44766_clean_otus_norm.txt", sep="\t", header=T, row=1)
soil_d <- vegdist(t(soil_norm), method="bray")
soil_meta <- read.table("data/soil/clean_map.txt", sep="\t", header=T)
soil_meta$ph_group <- cut(soil_meta$ph, breaks = c(0,4:8,14))
soil_cols <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black")[soil_meta$ph_group]
soil_bgs <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white")[soil_meta$ph_group]
soil_shps <- c(22, 25, 24, 21, 23, 23)[soil_meta$ph_group]
# Simulated dataset
sim_norm <- read.table("data/sim_gradient/fake_rel_abun_g1000_sd50_n50.txt", row=1, header=T, sep="\t")
sim_meta <- data.frame(SampleID=rownames(sim_n), noise=c(rep("n",50),rep("y",50)),
                       gradient=rep(round(seq(1,1000,length.out=50)),2))
sim_cols <- alpha(rep(viridis::viridis(50, alpha=0.8),2),0.6)
sim_d <- vegdist(sim_norm, nethod="bray")
# Cecum dataset
cecum_d <- readRDS("data/cecum/cecum_bray_dist.rds") # Bray-Curtis distances
cecum_meta <- read.table("data/cecum/cecum_meta.csv", sep=",", row=1, header=T)
## get rid of 7 samples removed for the paper
remove <- c("WPC2FHBD15B06C","WPC2FHBD18B07C","WPC2FHBD15B05C","WPC2FD15B04C",
            "WPC2FD15B05C","WPC2FD08B08C","WPC2FD01B10C")
remove_idx <- which(cecum_meta$SeqID %in% remove)
cecum_d <- as.dist(as.matrix(cecum_d)[-remove_idx,-remove_idx])
cecum_meta <- cecum_meta[-remove_idx,]
cecum_cols <- alpha(c("#fd2beb","#8e0013","#fa4d3f","#f99314","#00ae03","#57f89f","#0e2dfc","#6033f7"),0.6)[as.factor(meta_cecum$Collection)]
cecum_shp <- c(16, 17)[as.factor(cecum_meta$System)]


# Determine a reasonable degree:n ratio cutoff base donthese datasets
source('lib/lgd_source.r') # note that for this we have lg.dist return the ratio!

## Soil data - looks like 1:10 avgdeg:n ratio is reasonable cutoff here
try <- round(seq(max(soil_d), min(soil_d), length.out=20),3) # chose 0.585 because afterwords avg deg of graph is 2
par(mfrow=c(4,5), mar=c(3,3,4,1), mgp=c(1.5,0.5,0))
for (r in try) {
  print(paste0("radius: ",r))
  out <- lg.dist(soil_d, r)
  pc <- cmdscale(out[[1]], k=80, eig=T)
  pvar <- round((pc$eig/sum(pc$eig)) * 100, 1)
  plot(pc$points[,1], pc$points[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs, asp=1,
       main=paste0("radius ",r), xlab=paste0("PC 1 [",pvar[1],"%]"), ylab=paste0("PC 2 [",pvar[2],"%]"))
  legend("bottomleft", legend=paste0("ratio: ", out[[2]]), bg=alpha("white", 0.7),
         inset=0.2, x.intersp=0, y.intersp=0)
  print(" ")
}
## note: exported as 1100x800 to supplemental_figures

## Simulated data - looks like 1:10 ratio is ok here too
try <- round(seq(max(sim_d), min(sim_d), length.out=20),3)
par(mfrow=c(4,5), mar=c(3,3,4,1), mgp=c(1.5,0.5,0))
for (r in try) {
  print(paste0("radius: ",r))
  out <- lg.dist(sim_d, r)
  pc <- cmdscale(out[[1]], k=80, eig=T)
  pvar <- round((pc$eig/sum(pc$eig)) * 100, 1)
  plot(pc$points[,1], pc$points[,2], pch=20, col=sim_cols, main=paste0("radius ",r),
       xlab=paste0("PC 1 [",pvar[1],"%]"), ylab=paste0("PC 2 [",pvar[2],"%]"), asp=1)
  legend("bottomleft", legend=paste0("ratio: ", out[[2]]), bg=alpha("white", 0.7),
         inset=0.2, x.intersp=0.1, y.intersp=0)
  print(" ")
}
## note: exported as 1100x800 to supplemental_figures

## Cecum data - more sensitive, faster approach to the Y-shape
try <- round(seq(max(cecum_d), min(cecum_d), length.out=20),3)
par(mfrow=c(4,5), mar=c(3,3,4,1), mgp=c(1.5,0.5,0))
for (r in try) {
  print(paste0("radius: ",r))
  out <- lg.dist(cecum_d, r)
  pc <- cmdscale(out[[1]], k=80, eig=T)
  pvar <- round((pc$eig/sum(pc$eig)) * 100, 1)
  plot(pc$points[,1], pc$points[,2], pch=cecum_shp, col=cecum_cols, bg=cecum_cols,
       main=paste0("radius ",r), xlab=paste0("PC 1 [",pvar[1],"%]"), ylab=paste0("PC 2 [",pvar[2],"%]"), asp=1)
  legend("bottomleft", legend=paste0("ratio: ", out[[2]]), bg=alpha("white", 0.7),
         inset=0.2, x.intersp=0.1, y.intersp=0)
  print(" ")
}
## note: exported as 1100x800 to supplemental_figures



# Testing updates to the algorithm (can try multiple radii, optimization function)
source('lib/lgd_source.r')
## try running the default algorithm and see what is picked
soil_lgd <- lg.dist(soil_d)
## visually compare
soil_pc <- cmdscale(soil_d, k=10, eig=T)
pvar <- round((soil_pc$eig/sum(soil_pc$eig)) * 100, 1)
soil_pc_lgd <- cmdscale(soil_lgd, k=10, eig=T)
pvar_lgd <- round((soil_pc_lgd$eig/sum(soil_pc_lgd$eig)) * 100, 1)
par(mfrow=c(1,2))
plot(soil_pc$points[,1], soil_pc$points[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs,
     main=paste0("Original Soil"), xlab=paste0("PC 1 [",pvar[1],"%]"), ylab=paste0("PC 2 [",pvar[2],"%]"), asp=1)
plot(soil_pc_lgd$points[,1], soil_pc_lgd$points[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs,
     main=paste0("LGD-adjusted Soil"), xlab=paste0("PC 1 [",pvar_lgd[1],"%]"), ylab=paste0("PC 2 [",pvar_lgd[2],"%]"), asp=1)





