# Oral Preliminary Exam: Exploration
# Suzie Hoops
# December 2022

##### Set Up #####
library(vegan)
library(lsa)
source('lib/lgd_source.r')
set.seed(25)



##### HELPER FUNCTIONS #####
# In PCoA plots, get the percent variance
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}



##### Loading Data #####
# Soil dataset (88 soils)
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
soil_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
soil_bor <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black") # borders (col)
soil_fil <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white") # fills (bg)
soil_shp <- c(22, 25, 24, 21, 23, 23) # shapes
ord <- cut(meta_soil$ph, breaks = c(0,4:8,14)) # pH groups (<4, 4-5, 5-6, 6-7, 7-8, >8)
soil_bor <- soil_bor[ord]; soil_fil <- soil_fil[ord]; soil_shp <- soil_shp[ord]

##### Soil Data: Try every dist/diss & PCoA #####
dist_pcoa <- function (x, method="euclidean") {
  if (method == "cosine") {
    cos <- cosine(as.matrix(soil_n))
    d <- as.dist(1-cos)
  } else if (method == "aitchison") {
    d <- vegdist(t(x)+1, method=method)
  } else {
    d <- vegdist(t(x), method=method)
  }
  pc <- cmdscale(d, k=2, eig=T)
  plot(pc$points, col=soil_bor, bg=soil_fil, pch=soil_shp, cex=2, lwd=3,
       xlab=paste0("PC 1 [", calc.perc.var(pc$eig,1), "%]"),
       ylab=paste0("PC 1 [", calc.perc.var(pc$eig,2), "%]"), cex.lab=2)
  title(paste0("Dist: ", method))
}
dist_list <- c("bray","euclidean","manhattan","cosine","jaccard","aitchison","chisq")
par(mfrow=c(3,3))
for (dm in dist_list) {
  dist_pcoa(soil_n, method=dm)
}


##### Soil Data: Try every Ordination & "best" dist/diss #####
