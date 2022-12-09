# playing with bounds
library(vegan)

# matrix maximizing euclidean distance
mat <- matrix(c(0.5,0,0,0.5,0,1,0,0,0,0,1,0), ncol=4, byrow=T)
mat

# euclidean
mat_d <- vegdist(mat, method="euclidean")
mat_d
sqrt(2)    # relationship b/w 2 and 3
sqrt(1.5)  # relationship between 1 and 2,3

# aitchison
mat_ad1 <- vegdist(mat+1, method="aitchison") # note we cannot have 0s in CLR data
mat_ad1

mat_ad2 <- vegdist((mat*1000)+1, method="aitchison")
mat_ad2

mat_ad3 <- vegdist(mat+0.000000000001, method="aitchison")
mat_ad3


# sparse matrices
library(Matrix)
x <- round(runif(200, min=1, max=5000))
i <- round(runif(200, min=1, max=100))
j <- round(runif(200, min=1, max=2000))
my_spm <- sparseMatrix(i=i,j=j,x=x)
my_spm <- my_spm[rowSums(my_spm) > 1000,]
my_spm_n <- sweep(my_spm, 1, rowSums(my_spm), "/")
plot(density(vegdist(my_spm, method="euclidean")))
plot(density(vegdist(my_spm_n, method="euclidean")))


# limitation of LGD : multimodal density of dists?
source("lib/lgd_source.r")
## group1=first half, group2=second half
make_samp <- function(loc) {
  if(loc=="start"){
    samp <- round(c(runif(200, min=0, max=1000),runif(200, min=0, max=100)))
  }
  else{
    samp <- round(c(runif(200, min=0, max=100), runif(200, min=0, max=1000)))
  }
  return(samp)
}
## make dataset
mm_dat <- as.data.frame(matrix(make_samp("start"), nrow=1))
for(i in 1:49) {
  if ((i %% 2) == 0) {
    mm_dat <- rbind(mm_dat, make_samp("start"))
  } else{
    mm_dat <- rbind(mm_dat, make_samp("end"))
  }
}
mm_d <- vegdist(mm_dat, method="bray")
plot(density(mm_d))
mm_pc <- cmdscale(mm_d, k=2, eig=F)
plot(mm_pc, col=rep(c("blue","purple")), pch=16)
# lgd limit placed between peaks of dist distribution
mm_lgd <- lg.dist(mm_d, neighborhood.radius = 0.6)
mm_pc_lgd <- cmdscale(mm_lgd, k=2, eig=F)
plot(mm_pc_lgd, col=rep(c("blue","purple")), pch=17)
# more extremely low lgd placement
mm_lgd2 <- lg.dist(mm_d, neighborhood.radius = 0.3)
mm_pc_lgd2 <- cmdscale(mm_lgd2, k=2, eig=F)
plot(mm_pc_lgd2, col=rep(c("blue","purple")), pch=17)
