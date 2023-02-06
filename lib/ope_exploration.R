# Oral Preliminary Exam: Exploration & Images
# Suzie Hoops
# December 2022

##### Set Up #####
library(vegan)
library(lsa)
library(scales)
library(stringr)
library(Rtsne)
library(rdd)
source('lib/lgd_source.r')
set.seed(25)



##### HELPER FUNCTIONS #####
# Percent variance in PCoA plots
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}

# Ordination plot
ord_plot <- function (dat, mytitle="", mycols, mybor, myshp=16, flip1=F, flip2=F, eig=NULL) {
  # flip axes if needed
  if (flip1) dat[,1] <- dat[,1] * -1
  if (flip2) dat[,2] <- dat[,2] * -1
  # square plot output
  minlim <- min(dat[,1:2])
  maxlim <- max(dat[,1:2])
  # determine point size:
  if (nrow(dat) > 100) {sizepts <- 2}
  else {sizepts <- 3}
  # plot (add percentages if pcoa plot)
  if (is.null(eig)) {
    plot(dat, xlim=c(minlim, maxlim), ylim=c(minlim, maxlim),
         col=mycols, bg=mybor, pch=myshp, cex=sizepts, lwd=sizepts,
         xlab="Axis 1", ylab="Axis 2", cex.lab=1.5)
    title(mytitle, cex.main=2)
  } else {
    plot(dat, xlim=c(minlim, maxlim), ylim=c(minlim, maxlim),
         col=mycols, bg=mybor, pch=myshp, cex=sizepts, lwd=sizepts,
         xlab=paste0("PC1 [", calc.perc.var(eig,1), "%]"),
         ylab=paste0("PC2 [", calc.perc.var(eig,2), "%]"), cex.lab=1.5)
    title(mytitle, cex.main=2)
  }
}

# Plotting PCoA for given distance
dist_pcoa <- function (x, method="euclidean", mycols, mybor, myshp=16, flippc1=F, flippc2=F) {
  # distance calculation
  if (method == "cosine") {
    cos <- cosine(as.matrix(t(x)))
    d <- as.dist(1-cos)
  } else if (method == "aitchison") {
    d <- vegdist(x+0.00001, method=method)
  } else {
    d <- vegdist(x, method=method)
  }
  # calculate pc, flip axes if needed
  pc <- cmdscale(d, k=2, eig=T)
  # plot
  ord_plot(pc$points, mytitle=str_to_title(method), mycols, mybor, myshp, flippc1, flippc2, pc$eig)
}

# Smoothing: LGD
get_lgd <- function (d, grad=NULL, set_bestr=NULL) {
  ## Apply LGD with 25 values of r, get best r index & list of outputs
  ##  --> currently using best corr b/w as "best r"
  ##  --> should also try the biggest right skew
  print("Beginning LGD computation:")
  
  # initialize
  lgdlist <- list()
  rvals <- seq(min(d), max(d), length.out=25)
  bestidx <- 1   # index of best r
  bestval <- 0   # value for the best r
  
  # loop over r values
  for (i in 1:length(rvals)) {
    r <- rvals[i]
    print(paste0("  computing LGD for r = ", r))
    # computing LGD & save to list
    lg_d <- as.dist(lg.dist(d, neighborhood.radius = r))
    lgdlist <- append(lgdlist, list(lg_d))
    # determine best r value (default: correlation with PC dists, or use gradient)
    if (is.null(grad)) {
      pc_d <- vegdist(cmdscale(lg_d, k=3, eig=F), method="euclidean")
      curr_val <- r - round(cor(pc_d, lg_d, method="pearson"),4)
    } else {
      grad_d <- vegdist(grad, method="euclidean")
      curr_val <- round(cor(grad_d, lg_d, method="pearson"),4)
    }
    if (is.na(curr_val)) {curr_val <- 0}
    if (curr_val > bestval) {
      print(paste0("  --> new best r found, corr = ", curr_val))
      bestidx <- i
      bestval <- curr_val
    }
  }
  
  # overwrite best r if given as parameter:
  if (!is.null(set_bestr)) {
    found_bestidx <- bestidx
    bestidx <- which.min(abs(rvals - set_bestr))
    print(paste0("best r set to: ", round(rvals[bestidx],4), ", ",
                 round(rvals[bestidx] - rvals[found_bestidx],4)," away from best r found."))
  }
  
  # return: lgd list, best r index
  out <- list("lgdlist"=lgdlist, "bestidx"=bestidx, "rvals"=rvals)
  return(out)
}

# Smoothing: weighted results
get_wts <- function (rvals, bestidx, title="provided data") {
  # Determine weights centered on best r index provided, return weights
  wts <- kernelwts(1:length(rvals), center=bestidx, bw=4, kernel="gaussian")
  plot(rvals, wts, main=paste0("Weights: ", title), pch=16, cex=2.5, xlab="radius value", ylab="weight")
  return(wts)
}
get_wmeans <- function (lgdlist, bestidx, rvals, title="provided data") {
  # Determine new distances as weighted mean of other r distance outputs
  ## (1) get weights centered at best r index
  wts <- get_wts(rvals, bestidx, title)
  ## (2) per pairwise distance, get weighted mean
  means <- lgdlist[[1]]    # temporarily set means to first distance object
  for (i in 1:length(lgdlist[[1]])) {
    x <- sapply(lgdlist, function (l) { l[[i]] })
    means[i] <- weighted.mean(x, wts)
  }
  ## (3) return in a distance object
  return(as.dist(means))
}


##### Loading Data #####
# Simulated dataset (normal noise)
sim_n <- read.table("data/sim_gradient/fake_rel_abun_g1000_sd50_n50.txt", row=1, header=T, sep="\t")
meta_sim <- data.frame(SampleID=rownames(sim_n), noise=c(rep("n",50),rep("y",50)),
                       gradient=rep(round(seq(1,1000,length.out=50)),2))
sim_cols <- alpha(rep(viridis::viridis(50, alpha=0.8),2),0.6)

# Soil dataset (88 soils)
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
meta_soil$depth_mm <- as.numeric(gsub(".*-","",meta_soil$depth))
soil_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
## colors and shapes
soil_bor <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black")   # borders (col)
soil_fil <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white")   # fills (bg)
soil_shp <- c(22,25,24,21,23,23)   # shapes
ord <- cut(meta_soil$ph, breaks = c(0,4:8,14))   # pH groups (<4, 4-5, 5-6, 6-7, 7-8, >8)
soil_bor <- soil_bor[ord]; soil_fil <- soil_fil[ord]; soil_shp <- soil_shp[ord]

# Guerrero Negro dataset
meta_gn <- read.table("data/guerrero_negro/clean_map.txt", header=T, sep="\t")
rownames(meta_gn) <- meta_gn$SampleID
meta_gn$depth_mm <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-10","10-22","22-34")[as.factor(meta_gn$depth)]
gn_n <- read.table("data/guerrero_negro/47908_clean_otus_norm.txt", row=1, header=T, sep="\t")
## colors and shapes
depth_cols1 <- c("#F2182F", "#FF6A4F", "#FFA96A", "#FFDE95", "#FCF7C2",
                 "#D6F0F6", "#88D9E8", "#23AED1", "#0076B4")
ord_d <- as.factor(meta_gn$depth)
depth_cols <- depth_cols1[ord_d]

# Turkey Cecum dataset
cecum_d <- readRDS("data/cecum/cecum_bray_dist.rds") # Bray-Curtis distances
meta_cecum <- read.table("data/cecum/cecum_meta.csv", sep=",", row=1, header=T)
## get rid of 7 samples removed for the paper
remove <- c("WPC2FHBD15B06C","WPC2FHBD18B07C","WPC2FHBD15B05C","WPC2FD15B04C",
            "WPC2FD15B05C","WPC2FD08B08C","WPC2FD01B10C")
remove_idx <- which(meta_cecum$SeqID %in% remove)
cecum_d <- as.dist(as.matrix(cecum_d)[-remove_idx,-remove_idx])
meta_cecum <- meta_cecum[-remove_idx,]
## colors and shapes
cecum_col1 <- alpha(c("#fd2beb","#8e0013","#fa4d3f","#f99314","#00ae03","#57f89f","#0e2dfc","#6033f7"),0.6)
cecum_shp1 <- c(16, 17)
ord_c <- as.factor(meta_cecum$Collection)
ord_s <- as.factor(meta_cecum$System)
cecum_col = cecum_col1[ord_c]; cecum_shp = cecum_shp1[ord_s]

# MAGIC dataset
magic_n <- read.table("data/magic/mycleaned_magic_norm.txt", row=1, header=T, sep="\t")
meta_magic <- read.table("data/magic/clean_metafile.txt", header=T, sep="\t")
rownames(meta_magic) <- meta_magic$Sample_ID
meta_magic <- meta_magic[colnames(magic_n),]
meta_magic$Timeline_Weeks <- as.integer(meta_magic$Timeline_Weeks)
## refine dataset
meta_magic$age_grp <- cut(meta_magic$Timeline_Weeks, breaks = c(2,5,13,25,37,49,61,73,100)) # excluding first 2 weeks
remove_idx <- which(is.na(meta_magic$age_grp) | is.na(meta_magic$currentfeed_bf))           # refine to age set and breastfeeding info
meta_magic <- meta_magic[-remove_idx,]
magic_n <- magic_n[,-remove_idx]
## colors and shapes
magic_cols1 <- alpha(c("#ffb14e","#fa7b67","#cd34b5","#9d02d7","#0000ff","#009aa7","#40c557","#036003"),0.6)
magic_shp1 <- c(17,16)
ord_c <- as.factor(meta_magic$age_grp)
ord_s <- as.factor(meta_magic$currentfeed_bf)
magic_cols <- magic_cols1[ord_c]; magic_shp <- magic_shp1[ord_s];



##### Try every dist/diss with PCoA #####
# Distances
dist_list <- c("bray","kulczynski","euclidean","manhattan","cosine","jaccard","aitchison","chisq")

# Simulated data
par(mfrow=c(2,4),mgp=c(2.5, 1, 0))
for (dm in dist_list) {
  if (dm %in% c("bray","kulczynski","euclidean","manhattan")) {
    dist_pcoa(sim_n, method=dm, sim_cols, sim_cols, 16, flippc1=T)
  }
  else if (dm %in% c("aitchison")) {
    dist_pcoa(sim_n, method=dm, sim_cols, sim_cols, 16, flippc1=T, flippc2=T)
  }
  else if (dm %in% c("cosine")) {
    dist_pcoa(sim_n, method=dm, sim_cols, sim_cols, 16, flippc2=T)
  } else {
    dist_pcoa(sim_n, method=dm, sim_cols, sim_cols, 16)
  }
}
par(mfrow=c(1,1),mgp=c(3, 1, 0))
## note: ChiSq distance has greatest correlation (0.9934) and best results from visual inspection

# Soil data
par(mfrow=c(2,4),mgp=c(2.5, 1, 0))
for (dm in dist_list) {
  if (dm %in% c("bray","kulczynski","euclidean","manhattan","cosine","jaccard","aitchison")) {
    dist_pcoa(t(soil_n), method=dm, soil_bor, soil_fil, soil_shp, flippc2=T)
  } else {
    dist_pcoa(t(soil_n), method=dm, soil_bor, soil_fil, soil_shp)
  }
}
par(mfrow=c(1,1),mgp=c(3, 1, 0))
## note: most are correlated about the same, might want to pick jaccard (cor 0.924) b/c closest to paper

# Turkey cecum data - Can't do all distances because only have bray curtis distances object

# Guerrero Negro data
par(mfrow=c(2,4),mgp=c(2.5, 1, 0))
for (dm in dist_list) {
  if (dm %in% c("bray","manhattan","aitchison")) {
    dist_pcoa(t(gn_n), method=dm, depth_cols, depth_cols, 16, flippc2 = T)
  } else{
    dist_pcoa(t(gn_n), method=dm, depth_cols, depth_cols, 16)
  }
}
par(mfrow=c(1,1),mgp=c(3, 1, 0))


# MAGIC data: seems that Bray-Curtis is an ok place to start, but best is Aitchison
par(mfrow=c(2,4),mgp=c(2.5, 1, 0))
for (dm in dist_list) {
  dist_pcoa(t(magic_n), method=dm, magic_cols, magic_cols, 16)
}
par(mfrow=c(1,1),mgp=c(3, 1, 0))


##### Try every Ordination & "best" dist/diss #####
# Ordination approaches: PCoA, NMDS, t-SNE, constrained
# Simulated data (ChiSq distance)
## calculations
sim_chisq <- vegdist(sim_n, method="chisq")
sim_pc <- cmdscale(sim_chisq, k=2, eig=T)                           # PCoA
sim_nmds <- metaMDS(sim_chisq, k=2, try=20, trymax=50, maxit=500)   # NMDS
sim_tsne <- Rtsne(sim_chisq, is_distance=T, dims=3)                 # t-SNE
sim_rda <- dbrda(sim_n ~ gradient, distance="chisq", data=meta_sim) # constrained
(fit <- envfit(sim_rda, meta_sim[,2:3], perm = 999))
## plotting (note: dbRDa is slightly more complex)
par(mfrow=c(2,2),mgp=c(2.5, 1, 0))
ord_plot(sim_pc$points, "PCoA", sim_cols, sim_cols, 16, eig=sim_pc$eig)
ord_plot(sim_nmds$points, "NMDS", sim_cols, sim_cols, 16)
ord_plot(sim_tsne$Y, "t-SNE", sim_cols, sim_cols, 16)
plot(sim_rda, type="n", xlim=c(-1,1), ylim=c(-1,1), xlab="dbRDA 1", ylab="MDS 1", cex.lab=1.5)
points(sim_rda, col=sim_cols, bg=sim_cols, pch=16, cex=3)
arrows(0,0,1,0, lwd=2, length=0.1)
text(0.5, 0, pos=3, labels="gradient")
title("db-RDA", cex.main=2)
par(mfrow=c(1,1),mgp=c(3, 1, 0))


# Soil data (Jaccard dissimilarity)
soil_jacc <- vegdist(t(soil_n), method="jaccard")
soil_pc <- cmdscale(soil_jacc, k=2, eig=T)                              # PCoA
soil_nmds <- metaMDS(soil_jacc, k=2, try=20, trymax=50, maxit=500)      # NMDS
soil_tsne <- Rtsne(soil_jacc, is_distance=T, dims=3, perplexity=20)     # t-SNE
soil_rda <- dbrda(t(soil_n) ~ ph + latitude + annual_season_temp, distance="jaccard", data=meta_soil)      # constrained
(fit <- envfit(soil_rda, meta_soil[,c(45,31,43)], perm = 999))
par(mfrow=c(2,2),mgp=c(2.5, 1, 0))
ord_plot(soil_pc$points*-1, "PCoA", soil_bor, soil_fil, soil_shp)
ord_plot(soil_nmds$points, "NMDS", soil_bor, soil_fil, soil_shp)
ord_plot(soil_tsne$Y, "t-SNE", soil_bor, soil_fil, soil_shp)
plot(soil_rda, type="n", xlab="dbRDA 1", ylab="MDS 1", cex.lab=1.5)
points(soil_rda, col=soil_bor, bg=soil_fil, pch=soil_shp, cex=3)
text(soil_rda, dis="cn", cex=1)
title("db-RDA", cex.main=2)
par(mfrow=c(1,1),mgp=c(3, 1, 0))

# Turkey Cecum (Bray-Curtis dissimilarity)
cecum_pc <- cmdscale(cecum_d, k=2, eig=T)                               # PCoA
cecum_nmds <- metaMDS(cecum_d, k=2, try=20, trymax=50, maxit=500)       # NMDS
cecum_tsne <- Rtsne(cecum_d, is_distance=T, dims=3, perplexity=20)      # t-SNE
cecum_rda <- dbrda(cecum_d ~ System + Collection, data=meta_cecum)      # constrained
(fit <- envfit(cecum_rda, meta_cecum[,c(3,4)], perm = 999))
par(mfrow=c(2,2),mgp=c(2.5, 1, 0))
ord_plot(cecum_pc$points*-1, "PCoA", cecum_col, cecum_col, cecum_shp)
ord_plot(cecum_nmds$points, "NMDS", cecum_col, cecum_col, cecum_shp)
ord_plot(cecum_tsne$Y, "t-SNE", cecum_col, cecum_col, cecum_shp)
plot(cecum_rda, type="n", xlab="dbRDA 1", ylab="MDS 1", cex.lab=1.5)
points(cecum_rda, col=cecum_col, bg=cecum_col, pch=cecum_shp, cex=3)
text(cecum_rda, dis="cn", cex=1)
title("db-RDA", cex.main=2)
par(mfrow=c(1,1),mgp=c(3, 1, 0))

# Guerrero Negro (Jaccard dissimilarity)
gn_d <- vegdist(t(gn_n), method="jaccard")
gn_pc <- cmdscale(gn_d, k=2, eig=T)                               # PCoA
gn_nmds <- metaMDS(gn_d, k=2, try=20, trymax=50, maxit=500)       # NMDS
gn_tsne <- Rtsne(gn_d, is_distance=T, dims=3, perplexity=5)       # t-SNE
gn_rda <- dbrda(gn_d ~ end_depth + layer, data=meta_gn)           # constrained
(fit <- envfit(gn_rda, meta_gn[,c(30,40)], perm = 999))
par(mfrow=c(2,2),mgp=c(2.5, 1, 0))
ord_plot(gn_pc$points*-1, "PCoA", depth_cols, depth_cols, 16, flip1 = T, flip2=T)
ord_plot(gn_nmds$points, "NMDS", depth_cols, depth_cols, 16)
ord_plot(gn_tsne$Y, "t-SNE", depth_cols, depth_cols, 16)
plot(gn_rda, type="n", xlab="dbRDA 1", ylab="MDS 1", cex.lab=1.5)
points(gn_rda, col=depth_cols, bg=depth_cols, pch=16, cex=3)
text(gn_rda, dis="cn", cex=1)
title("db-RDA", cex.main=2)
par(mfrow=c(1,1),mgp=c(3, 1, 0))

# MAGIC (Aitchison distance)
magic_d <- vegdist(t(magic_n)+0.000001, method="aitchison")
magic_pc <- cmdscale(magic_d, k=2, eig=T)                               # PCoA
magic_nmds <- metaMDS(magic_d, k=2, try=20, trymax=50, maxit=500)       # NMDS
magic_tsne <- Rtsne(magic_d, is_distance=T, dims=3, perplexity=30)      # t-SNE
magic_rda <- dbrda(magic_d ~ Timeline_Weeks + currentfeed_bf, data=meta_magic) # constrained
(fit <- envfit(magic_rda, meta_magic[,c(6,650)], perm = 999))
par(mfrow=c(2,2),mgp=c(2.5, 1, 0))
ord_plot(magic_pc$points*-1, "PCoA", magic_cols, magic_cols, magic_shp, flip1=T)
ord_plot(magic_nmds$points, "NMDS", magic_cols, magic_cols, magic_shp)
ord_plot(magic_tsne$Y, "t-SNE", magic_cols, magic_cols, magic_shp)
plot(magic_rda, type="n", xlab="dbRDA 1", ylab="MDS 1", cex.lab=1.5)
points(magic_rda, col=magic_cols, bg=magic_cols, pch=magic_shp, cex=3)
text(magic_rda, dis="cn", cex=1)
title("db-RDA", cex.main=2)
par(mfrow=c(1,1),mgp=c(3, 1, 0))



##### Uneven Gradient #####
# Select uneven samples
select_uneven <- c(seq(1,10,2),15,20,25,30,35,seq(40,50,2))
select_uneven <- c(select_uneven, select_uneven+50)
uneven <- sim_n[select_uneven,]
uneven_d <- vegdist(uneven, method="bray")
# Visualize
uneven_pc <- cmdscale(uneven_d, k=2, eig=T)
plot(uneven_pc$points, col=sim_cols[select_uneven], pch=16, cex=3,
     xlab=paste0("PC 1 [",calc.perc.var(uneven_pc$eig, 1),"%]"),
     ylab=paste0("PC 1 [",calc.perc.var(uneven_pc$eig, 2),"%]"),
     main="Uneven Sampling")
# try LGD w/ smoothing
uneven_out <- get_lgd(uneven_d, meta_sim$gradient[select_uneven]) # best r found: 0.8 (maximum)
uneven_wtd <- get_wmeans(uneven_out$lgdlist, uneven_out$bestidx, uneven_out$rvals, "uneven gradient")
uneven_lgd_pc <- cmdscale(uneven_wtd, k=2, eig=T)
lims=c(min(uneven_lgd_pc$points), max(uneven_lgd_pc$points))
plot(uneven_lgd_pc$points, col=sim_cols[select_uneven], pch=16, cex=3,
     xlab=paste0("PC 1 [",calc.perc.var(uneven_lgd_pc$eig, 1),"%]"),
     ylab=paste0("PC 1 [",calc.perc.var(uneven_lgd_pc$eig, 2),"%]"),
     main="Uneven Sampling w/ LGD smoothing", xlim=lims, ylim=lims)
# try LGD w/o smoothing
uneven_lgd <- lg.dist(uneven_d, neighborhood.radius=0.8)
uneven_lgd_pc <- cmdscale(uneven_lgd, k=2, eig=T)
lims=c(min(uneven_lgd_pc$points), max(uneven_lgd_pc$points))
plot(uneven_lgd_pc$points, col=sim_cols[select_uneven], pch=16, cex=3,
     xlab=paste0("PC 1 [",calc.perc.var(uneven_lgd_pc$eig, 1),"%]"),
     ylab=paste0("PC 1 [",calc.perc.var(uneven_lgd_pc$eig, 2),"%]"),
     main="Uneven Sampling w/ LGD (r=0.8)", xlim=lims, ylim=lims)



##### Swiss Roll #####
# Loading Swiss Roll dataset (created with sklearn)
swiss_data <- read.delim("data/swiss_roll/swiss_roll.txt", sep=" ", header=F)
swiss_color <- c(read.delim("data/swiss_roll/swiss_roll_colors.txt", sep=" ", header=F)[,1])
swiss_euc <- vegdist(swiss_data, method="euclidean")
# Traditional PCA
swiss_pc <- cmdscale(swiss_euc, k=2, eig=F)
plot(swiss_pc, col=swiss_color, pch=16, cex=2, xlab="PC 1", ylab="PC 2", main="Swiss Roll PCA")
# Manifold Learning
swiss_tsne <- Rtsne(swiss_euc, dims=2, perplexity=45, initial_dims=3)
plot(swiss_tsne$Y, col=swiss_color, pch=16, cex=2, xlab="Axis 1", ylab="Axis 2", main="Swiss Roll t-SNE")
swiss_isomap <- isomap(swiss_euc, ndim=2, k=5)
plot(swiss_isomap$points, col=swiss_color, pch=16, cex=2, xlab="Axis 1", ylab="Axis 2", main="Swiss Roll Isomap")
# PCA w/ LGD
plot(density(swiss_euc), lwd=3, xlab="Distances", ylab="Density", main="Density of Swiss Dists")
swiss_lgd <- cmdscale(lg.dist(swiss_euc, neighborhood.radius=6.2), k=2, eig=F)
plot(swiss_lgd, col=swiss_color, pch=16, cex=2, xlab="PC 1", ylab="PC 2", main="Swiss Roll LGD")
##### Using LGD #####
par(mgp=c(2.5,1,0))
# simulated data
sim_out <- get_lgd(sim_chisq, meta_sim$gradient) # best r: 2.078
sim_wtd <- get_wmeans(sim_out$lgdlist, sim_out$bestidx, sim_out$rvals, "simulated data")
sim_pc <- cmdscale(sim_chisq, k=2, eig=T)
sim_ord1 <- cmdscale(lg.dist(sim_chisq, neighborhood.radius=1.65), k=2, eig=T)
sim_ord <- cmdscale(sim_wtd, k=2, eig=T)
ord_plot(sim_pc$points, "Simulated", sim_cols, sim_cols, 16, eig=sim_pc$eig)
ord_plot(sim_ord1$points, "Simulated (one radius)", sim_cols, sim_cols, 16, eig=sim_ord1$eig)
ord_plot(sim_ord$points, "Simulated (smooth LGD)", sim_cols, sim_cols, 16, eig=sim_ord$eig)
# soil data
soil_out <- get_lgd(soil_jacc, meta_soil$ph) # best r: 0.9804
soil_wtd <- get_wmeans(soil_out$lgdlist, soil_out$bestidx, soil_out$rvals, "soil data")
soil_pc <- cmdscale(soil_jacc, k=2, eig=T); soil_pc$points[,2] <- soil_pc$points[,2] * -1;
soil_ord1 <- cmdscale(lg.dist(soil_jacc, neighborhood.radius=0.9411), k=2, eig=T); soil_ord1$points[,2] <- soil_ord1$points[,2] * -1;
soil_ord <- cmdscale(soil_wtd, k=2, eig=T)
ord_plot(soil_pc$points, "Soil", soil_bor, soil_fil, soil_shp, eig=soil_pc$eig)
ord_plot(soil_ord1$points, "Soil (one radius)", soil_bor, soil_fil, soil_shp, eig=soil_ord1$eig)
ord_plot(soil_ord$points, "Soil (smooth LGD)", soil_bor, soil_fil, soil_shp, eig=soil_ord$eig)
cor.test(soil_ord1$points[,2], meta_soil$annual_season_temp) # cor: 0.597, p<0.0001
cor.test(soil_ord1$points[,2], meta_soil$latitude) # cor: -0.544, p < 0.0001
# guerrero negro data
gn_out <- get_lgd(vegdist(t(gn_n), method="jaccard"), meta_gn$end_depth, set_bestr=0.74) # best r found: 0.6403
gn_wtd <- get_wmeans(gn_out$lgdlist, gn_out$bestidx, gn_out$rvals, "guerrero negro")
gn_pc <- cmdscale(vegdist(t(gn_n), method="jaccard"), k=2, eig=T)
gn_ord1 <- cmdscale(lg.dist(vegdist(t(gn_n), method="jaccard"), neighborhood.radius=0.74), k=2, eig=T)
gn_ord <- cmdscale(gn_wtd, k=2, eig=T)
ord_plot(gn_pc$points, "Guerrero Negro", depth_cols, depth_cols, 16, eig=gn_pc$eig)
ord_plot(gn_ord1$points, "Guerrero Negro (one radius)", depth_cols, depth_cols, 16, eig=gn_ord1$eig)
ord_plot(gn_ord$points, "Guerrero Negro (smooth LGD)", depth_cols, depth_cols, 16, eig=gn_ord$eig)
# turkey cecum data
cecum_out <- get_lgd(cecum_d) # best r: 0.4931
cecum_wtd <- get_wmeans(cecum_out$lgdlist, cecum_out$bestidx, cecum_out$rvals, "turkey cecum")
cecum_pc <- cmdscale(cecum_d, k=2, eig=T); cecum_pc$points[,2] <- cecum_pc$points[,2] * -1
cecum_ord1 <- cmdscale(lg.dist(cecum_d, neighborhood.radius=0.93), k=2, eig=T)
cecum_ord <- cmdscale(cecum_wtd, k=2, eig=T); cecum_ord$points[,2] <- cecum_ord$points[,2] * -1;
ord_plot(cecum_pc$points, "Turkey Cecum", cecum_col, cecum_col, cecum_shp, eig=cecum_pc$eig)
ord_plot(cecum_ord1$points, "Turkey Cecum (one radius)", cecum_col, cecum_col, cecum_shp, eig=cecum_ord1$eig)
ord_plot(cecum_ord$points, "Turkey Cecum (smooth LGD)", cecum_col, cecum_col, cecum_shp, eig=cecum_ord$eig)
# MAGIC data (Bray-Curtis)
magic_bc <- vegdist(t(magic_n), method="bray")
magic_out <- get_lgd(magic_bc)
magic_wtd <- get_wmeans(magic_out$lgdlist, magic_out$bestidx, magic_out$rvals, "MAGIC")
magic_pc <- cmdscale(magic_bc, k=2, eig=T); magic_pc$points[,1] <- magic_pc$points[,1] * -1;
magic_ord1 <- cmdscale(lg.dist(magic_bc, neighborhood.radius=0.93), k=2, eig=T); magic_ord1$points[,1] <- magic_ord1$points[,1] * -1;
magic_ord <- cmdscale(magic_wtd, k=2, eig=T)
ord_plot(magic_pc$points, "MAGIC", magic_cols, magic_cols, magic_shp, eig=magic_pc$eig)
ord_plot(magic_ord1$points, "MAGIC (one radius)", magic_cols, magic_cols, magic_shp, eig=magic_ord1$eig)
ord_plot(magic_ord$points, "MAGIC (smooth LGD)", magic_cols, magic_cols, magic_shp, eig=magic_ord$eig)
# MAGIC data (Aitchison)
magic_a <- vegdist(t(magic_n)+0.00001, method="aitchison")
magic_out <- get_lgd(magic_a)
magic_wtd <- get_wmeans(magic_out$lgdlist, magic_out$bestidx, magic_out$rvals, "MAGIC")
magic_pc <- cmdscale(magic_a, k=2, eig=T); magic_pc$points[,1] <- magic_pc$points[,1] * -1;
magic_ord1 <- cmdscale(lg.dist(magic_a, neighborhood.radius=100), k=2, eig=T)
magic_ord <- cmdscale(magic_wtd, k=2, eig=T); magic_ord$points[,2] <- magic_ord$points[,2] * -1;
ord_plot(magic_pc$points, "MAGIC", magic_cols, magic_cols, magic_shp, eig=magic_pc$eig)
ord_plot(magic_ord1$points, "MAGIC (one radius)", magic_cols, magic_cols, magic_shp, eig=magic_ord1$eig)
ord_plot(magic_ord$points, "MAGIC (smooth LGD)", magic_cols, magic_cols, magic_shp, eig=magic_ord$eig)



##### Sensitivity Analysis cont. #####
# Lengthen gradient until distances become oversaturated and arch appears

# Simulated dataset generation (no noise for now)
## settings: gradient length maximum 1000, std dev 50, 50 samples
m <- 1000 # max gradient length
sd <- 50  # standard dev
n <- 50   # num samples
## create coenoclines
xvals <- seq(-(4*sd), m+(4*sd), by=sd/4)[-1]                     # ensures 0 to m in figure is flat
x <- sapply(xvals, function (x) {dnorm(seq(0,m,1), x, sd)*100})  # create null distributions
x[abs(x) < 0.001] <- 0                                           # remove otus entirely at tails
x <- sweep(x, 1, rowSums(x), "/")                                # normalize rows (i.e. possible sample combos)
x <- x[,-which(colSums(x) == 0)]                                 # remove buffer otus
## sample from the coenoclines
create_sim_data <- function (g) {
  # evenly sample gradient up to g provided
  inds <- round(seq(0, g, length.out=(n+2)))[-c(1,n+2)]
  # format table to relative abundances
  dat <- as.data.frame(x[inds,])
  dat <- dat[,-which(colSums(dat) == 0)]
  rownames(dat) <- paste0("sample.", rownames(dat))
  colnames(dat) <- paste0("otu.", 1:ncol(dat))
  # get colors
  dat_colors <- viridis::viridis(m)[inds]
  # return completed table and sample colors
  return(list(data=dat, colors=dat_colors))
}

# Plots: density of distances & pcoa
library(manipulate)
par(mfrow=c(1,2),mgp=c(2.5,1,0))
intract_plots <- function (g) {
  # calculations
  mysim <- create_sim_data(g)
  d <- vegdist(mysim$data, method="bray")
  pc <- cmdscale(d, k=2)
  minax <- min(pc); maxax <- max(pc);
  # plots
  plot(density(d), col="blue", lwd=4, xlab="Distance", ylab="Density", main="")
  plot(pc, col=alpha(mysim$colors,0.9), xlim=c(minax, maxax), ylim=c(minax, maxax),
       pch=16, cex=2.5, xlab="PC 1", ylab="PC 2", main="")
}
manipulate(intract_plots(g), g=slider(n,m,step=50))
par(mfrow=c(1,1),mgp=c(3,1,0))



##### SENS. ANALYSIS SHINY APP #####
# Try again with Shiny to make it cleaner.
library(shiny)
# Define UI
ui <- fluidPage(
  # Title
  titlePanel("Sliding length of simulated gradient."),
  # Input: Gradient Length
  sliderInput("g", "Gradient length:",
              min = 25, max = 1000,
              value = 25, step=25),
  # Display outputs
  plotOutput("plots")
)
# Define server
server <- function(input, output) {
  # Plot density & PCoA
  output$plots <- renderPlot({
    # calculations
    mysim <- create_sim_data(input$g)
    d <- vegdist(mysim$data, method="euclidean")
    pc <- cmdscale(d, k=2)
    minax <- min(pc); maxax <- max(pc);
    # plots
    op <- par(mfrow=c(1,2))
    plot(density(d), col="black", lwd=4, xlab="Distance", ylab="Density", main="Density of Distances")
    plot(pc, col=alpha(mysim$colors,0.9), xlim=c(minax, maxax), ylim=c(minax, maxax),
         pch=16, cex=2.5, xlab="PC 1", ylab="PC 2", main="PCoA")
    par(op)
  })
}

shinyApp(ui, server)





##### Scratch Stuff #####

# Mini simulated data example of smoothing
# par(mfrow=c(1,1),mgp=c(1, 1, 0))
# ## original
# mini_sim <-read.table("data/old_gradient/fake_rel_abun_long_n100.txt", row=1, header=T, sep="\t")
# mini_sim <- mini_sim[c(1,10,20,30,40,50,60,70,80,90,100),]
# mini_d <- vegdist(mini_sim, method="euclidean")
# mini_pc <- cmdscale(mini_d, k=2, eig=T)
# mini_cols <- viridis::viridis(11, alpha=0.8)
# mini_pc$points[,2] <- mini_pc$points[,2] * -1
# ## original plot
# plot(mini_pc$points, xaxt='n', yaxt='n',
#      xlim=range(mini_pc$points)+c(-0.02,0.02), ylim=range(mini_pc$points)+c(-0.01,0.02),
#      col=mini_cols, pch=16, cex=6,
#      xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# ## low r value plot
# mini_lgd1 <- lg.dist(mini_d, neighborhood.radius = 0.10)
# mini_lgd1_pc <- cmdscale(mini_lgd1, k=2, eig=T)
# plot(mini_lgd1_pc$points, xaxt='n', yaxt='n',
#      xlim=range(mini_lgd1_pc$points)+c(-0.02,0.02), ylim=range(mini_lgd1_pc$points)+c(-0.01,0.02),
#      col=mini_cols, pch=16, cex=6,
#      xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# ## medium r value plot
# mini_lgd2 <- lg.dist(mini_d, neighborhood.radius = 0.20)
# mini_lgd2_pc <- cmdscale(mini_lgd2, k=2, eig=T)
# plot(mini_lgd2_pc$points, xaxt='n', yaxt='n',
#      xlim=range(mini_lgd2_pc$points)+c(-0.02,0.02), ylim=range(mini_lgd2_pc$points)+c(-0.01,0.02),
#      col=mini_cols, pch=16, cex=6,
#      xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# ## high r value plot (use the above original one)
# ## smoothed result
# mini_out <- get_lgd(mini_d, 1:11) # best r: 0.1635 w/ corr 0.9996
# mini_wtd <- get_wmeans(mini_out$lgdlist, mini_out$bestidx, mini_out$rvals, "small simulated data")
# mini_smooth_pc <- cmdscale(mini_wtd, k=2, eig=T)
# plot(mini_smooth_pc$points, xaxt='n', yaxt='n',
#      xlim=range(mini_smooth_pc$points)+c(-0.02,0.02), ylim=range(mini_smooth_pc$points)+c(-0.01,0.02),
#      col=mini_cols, pch=16, cex=6,
#      xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# ## reset
# par(mfrow=c(1,1),mgp=c(3, 1, 0))