# Playing around with weighting results of LGD
# Effectively smoothing results
library(rdd)
source("lib/lgd_source.r")
set.seed(25)

##### Helper Functions #####
# Applying LGD to multiple r values
get_lgd <- function (d, grad) {
  ## Apply LGD with 25 values of r, get best r index & list of outputs
  ##  --> currently using best association with real gradient as "best r"
  ##  --> should also try the biggest right skew
  print("Beginning LGD computation:")
  
  # initialize
  lgdlist <- list()
  rvals <- seq(min(d), max(d), length.out=25)
  bestidx <- 1   # index of best r
  bestval <- 0   # value for the best r 
  grad_dists <- vegdist(grad, method="euclidean")
  
  # loop over r values
  for (i in 1:length(rvals)) {
    r <- rvals[i]
    print(paste0("  computing LGD for r = ", r))
    # computing LGD & save to list
    lg_d <- as.dist(lg.dist(d, neighborhood.radius = r))
    lgdlist <- append(lgdlist, list(lg_d))
    # determine best r value
    curr_val <- round(cor(grad_dists, lg_d, method="pearson"),4)
    if (is.na(curr_val)) {curr_val <- 0}
    if (curr_val >= bestval) {
      print(paste0("  --> new best r found, corr w/ gradient = ", curr_val))
      bestidx <- i
      bestval <- curr_val
    }
  }
  
  # return: lgd list, best r index
  out <- list("lgdlist"=lgdlist, "bestidx"=bestidx, "rvals"=rvals)
  return(out)
}

# In PCoA plots, get the percent variance
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}

# Plotting PCoA
plot_pcoa <- function (d, mycols, mybor, myshp=16, mytitle) {
  # calculate pcoa
  pc <- cmdscale(d, k=2, eig=T)
  # plot
  plot(pc$points,
       col=mycols, bg=mybor, pch=myshp, cex=2, lwd=3,
       xlab=paste0("PC1 [", calc.perc.var(pc$eig,1), "%]"),
       ylab=paste0("PC2 [", calc.perc.var(pc$eig,2), "%]"), cex.lab=2)
  title(mytitle, adj=0, cex.main=2)
}


##### Set Up & Data Loading #####
# simulated data - Chi-Sq distances
sim_n <- read.table("data/sim_gradient/fake_rel_abun_g1000_sd50_n50.txt", row=1, header=T, sep="\t")
meta_sim <- data.frame(SampleID=rownames(sim_n), noise=c(rep("n",50),rep("y",50)),
                       gradient=rep(round(seq(1,1000,length.out=50)),2))
sim_cols <- rep(viridis::viridis(50, alpha=0.8),2)
sim_d <- vegdist(sim_n, method="euclidean")

# soil data - Jaccard distances
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
soil_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
soil_d <- vegdist(t(soil_n), method="jaccard")
soil_col1 <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black") # colors (borders)
soil_bg1 <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white") # fills (background)
soil_shp1 <- c(22, 25, 24, 21, 23, 23) # shapes
ord <- cut(meta_soil$ph, breaks = c(0,4:8,14)) # pH groups (<4, 4-5, 5-6, 6-7, 7-8, >8)
soil_cols <- soil_col1[ord]; soil_bors <- soil_bg1[ord]; soil_shps <- soil_shp1[ord];


# LGD applied to simulated data (see helper function for current implementation)
sim_out <- get_lgd(sim_d, meta_sim$gradient) # best r: 0.187 w/ corr 0.9993
sim_lgdlist <- sim_out$lgdlist
sim_bestidx <- 2 # sim_out$bestidx
sim_rvals <- sim_out$rvals

# LGD applied to soil data (see helper function for current implementation)
soil_out <- get_lgd(soil_d, meta_soil$ph) # best r: 0.9411 w/ corr 0.7212
soil_lgdlist <- soil_out$lgdlist
soil_bestidx <- soil_out$bestidx
soil_rvals <- soil_out$rvals


##### Gaussian weight #####
# Function(s) for applying a gaussian weight to get output
get_wts <- function (rvals, bestidx) {
  # Determine weights centered on best r index provided, return weights
  wts <- kernelwts(1:length(rvals), center=bestidx, bw=4, kernel="gaussian")
  plot(rvals, wts, main="Weights", pch=16, cex=3)
  return(wts)
}
get_wmeans <- function (lgdlist, bestidx, rvals) {
  # Determine new distances as weighted mean of other r distance outputs
  ## (1) get weights centered at best r index
  wts <- get_wts(rvals, bestidx)
  ## (2) per pairwise distance, get weighted mean
  means <- lgdlist[[1]]    # temporarily set means to first distance object
  for (i in 1:length(lgdlist[[1]])) {
    x <- sapply(lgdlist, function (l) { l[[i]] })
    means[i] <- weighted.mean(x, wts)
  }
  ## (3) return in a distance object
  return(as.dist(means))
}

# Calculate weighted means of all pairwise distances
wtd_sim <- get_wmeans(sim_lgdlist, sim_bestidx, sim_rvals)
wtd_soil <- get_wmeans(soil_lgdlist, soil_bestidx, soil_rvals)



##### Plots comparing outputs #####
# Simulated Data
## 3 plots : no lgd, best lgd, smoothed lgd
par(mfrow=c(1,3))
plot_pcoa(sim_d, sim_cols, sim_cols, 16, "Simulated Data")
plot_pcoa(as.dist(lg.dist(sim_d, neighborhood.radius=sim_rvals[sim_bestidx])),
          sim_cols, sim_cols, 16, "LGD (one radius)")
plot_pcoa(wtd_sim, sim_cols, sim_cols, 16, "Smoothed LGD")
par(mfrow=c(1,1))

# Soil Data
## 3 plots : no lgd, best lgd, smoothed lgd
par(mfrow=c(1,3))
plot_pcoa(soil_d, soil_cols, soil_bors, soil_shps, "Simulated Data")
plot_pcoa(as.dist(lg.dist(soil_d, neighborhood.radius=soil_rvals[soil_bestidx])),
          soil_cols, soil_bors, soil_shps, "LGD (one radius)")
plot_pcoa(wtd_soil, soil_cols, soil_bors, soil_shps, "Smoothed LGD")
par(mfrow=c(1,1))
# running into a problem here! For some reason the soil data set has infinte values
#.  discovered one of the problems was a bad sample (BB1) which only had one sequence count,
#.  but now maybe HI1 is having a problem too and IDK why... shouldn't get infinite
#.  values from this...



##### Testing / Playing Around #####
# values to get a weighted mean from
myvals <- c(4,7,5,10,5,8,2,5,9,10,11,12,11,11,10,13,11,12,15,17) # outputs per r
# Create weights centered at the "best" index
mybestidx <- 14
mywts <- kernelwts(1:length(myvals), center=mybestidx, bw=4, kernel="gaussian")
plot(mywts, main="Weights")
plot(myvals, main="Values")
mymn <- weighted.mean(myvals, mywts)
abline(h=mymn, col="blue", lwd=3, lty="dashed")



# For overwriting soil data
tmp_m <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
tmp_m <- tmp_m[-which(tmp_m$SampleID %in% c("BB1")),]
write.table(tmp_m,"data/soil/clean_map.txt", sep="\t", col.names=T, row.names=F)

tmp_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
tmp_n <- tmp_n[,tmp_m$SampleID]
write.table(tmp_n, "data/soil/44766_clean_otus_norm.txt", sep="\t", col.names=T, row.names=T)

tmp_c <- read.table("data/soil/44766_clean_otus.txt", row=1, header=T, sep="\t")
tmp_c <- tmp_c[,tmp_m$SampleID]
write.table(tmp_c, "data/soil/44766_clean_otus.txt", sep="\t", col.names=T, row.names=T)


