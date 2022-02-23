# LLE_simulation.r
# Locally Linear Embedding on Simulation Data
# Knights Lab - University of Minnesota
# June 2019
# usage : LLE_simulation.r -i input_table.txt

##### Set Up #####
#library(optparse)
library(lle)


##### Parse Command Line #####
# option_list <- list(make_option(c("-i", "--input_table"), type="character",
#                                 help="Path to input file. Expects otu table in tab-delimited txt file format.") )
# opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun.txt", sep="\t", header=T)


##### Helpful Functions #####
colorfunc <- colorRampPalette(c("white", "navy")) # scale_color_manual(values=colorfunc(nrow(dat)))


##### Run Locally Linear Embedding #####
# Calculate Ideal Number of Neighbors
start_time <- Sys.time()
k <- calc_k(dat, m=1, kmin=1, kmax=25, T, T, 4)
end_time <- Sys.time()
print(paste0("time elapsed in k calculation: ", end_time-start_time, " seconds"))
k <- k$k[which.min(k$rho)]
# Run Local Linear Embedding Algorithm to create new coordinate field
start_time <- Sys.time()
lle_out <- lle(dat, m=2, k=k, reg=1, ss=F, id=T, v=0.9) # expects samples in rows
end_time <- Sys.time()
print(paste0("time elapsed in locally linear embedding: ", end_time-start_time, " seconds"))
lle_dims <- as.data.frame(lle_out$Y); colnames(lle_dims) <- c("dim1", "dim2");
ggplot(data.frame(X=lle_dims$dim1,Y=lle_dims$dim2, color=rep(1:200,2), noise=c(rep("no",200), rep("yes",200))),
       aes(x=X, y=Y, color=factor(color), shape=factor(noise), size=factor(noise))) +
  geom_point(alpha=0.6) +
  scale_color_viridis_d(nrow(lle_dims)/2) +
  scale_shape_manual(values=c(1, 18)) +
  scale_size_manual(values=c(15, 10)) +
  labs(title="LLE embedded data",
       x="first dimension",
       y="second dimension",
       shape="noise") +
  guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
  theme_classic() + NULL
plot( lle_out$id, main="intrinsic dimension", type="l", xlab=expression(x[i]), ylab="id", lwd=2 )


