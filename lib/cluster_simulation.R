# Generating a group-wise clusters simulation
# Suzie Hoops
# Last Updated: March 2023

library(stats)
library(plot.matrix)
library(vegan)
library(scales)
library(gtools)
set.seed(125)


# playing around with the cluster
# what do real samples look like? Let's consider the Soil dataset
soil <- read.delim("data/soil/44766_clean_otus_norm.txt", sep="\t", row=1, header=T) # loading soil dataset
soil <- soil[,runif(10, min=1, max=ncol(soil))]
soil <- soil[rowSums(soil) > 0,]
plot(density(soil[,1]), main="Soil Samples", type="l", lwd=1, col=alpha("blue",0.4), ylim=c(0,800))
for (c in 2:10) {
  lines(density(soil[,c]), lwd=1, col=alpha("blue",0.4))
}
legend(x=0.01, y=450, legend=c("Sample Distribution"), lwd=1, col="blue")
# can we emulate this with a beta distribution?
len_x <- seq(0,1,length.out=50)
prior <- dbeta(len_x, 1, 20)
plot(len_x, prior, type="l", lwd=3, main="Beta Distribution & Dirichlet Samples",
     xlab="Proportion of Sample", ylab="Frequency") # curve of prior distribution used
mysample <- rbeta(len_x, 1, 20) # sample the same curve
lines(density(mysample), lty="dashed", lwd=3, col="purple") # distribution of sample from beta
# what about other samples if we use this as input to dirichlet?
conf <- 10
mydirichlet <- rdirichlet(20, mysample*conf) # 20 dirichlet samples
mydirichlet[mydirichlet < .Machine$double.eps] <- 0
rowSums(mydirichlet)
for (row in 1:20) {
  lines(density(mydirichlet[row,],n=16), pch=16, cex=2, col=alpha("purple",0.2))
}
legend(x=0.4, y=12, legend=c("Beta(1,20)","Random Beta Sample","Dirichlet Samples"),
       lty=c("solid","dashed","solid"), lwd=c(3,3,1), col=c("black","purple","purple"))
# make a sample and visualize it
myset <- rbind(mysample/sum(mysample), mydirichlet)
heatmap(myset, col=topo.colors(254), Rowv=NA)
legend(x="right", legend=c("low","  |","  |","  |","high"), fill=topo.colors(5))



##### Function to create clusters dataset #####
# Create nclust clusters from nsamp samples
## note: balanced group sizes, but could add another parameter for inbalanced groups
create_clusters <- function (minsamp=90, nclust=3, n_otu=500, beta_ab=c(1,20)) {
  # What to do if overlap = T vs overlap = F??
  # Initialize variables needed
  nsamp_per <- ceiling(minsamp / nclust)    # number of samples per cluster
  len_x <- seq(0, 1, length.out=n_otu)      # x values for the beta function
  conf <- 10                                # confidence factor for dirichlet, currently built-in
  matrix <- matrix(NA, nrow=0, ncol=n_otu)  # empty matrix to fill with samples
  # Use beta distribution as prior to dirichlet of samples
  for (c in 1:nclust) {
    prior <- rbeta(len_x, beta_ab[1], beta_ab[2])
    dir_samples <- rdirichlet(nsamp_per, prior * conf)
    dir_samples[dir_samples < .Machine$double.eps] <- 0 # below machine epsilon set to 0
    matrix <- rbind(matrix, dir_samples)
  }
  # Clean up matrix and return
  matrix <- matrix[rowSums(matrix) > 0,]  # remove empty rows (should not happen, but safety check)
  matrix <- matrix[,colSums(matrix) > 0]  # remove empty columns
  return(matrix)
}
# Plotting output
mymat <- create_clusters()
heatmap(mymat, Rowv=NA, col=c("white",topo.colors(49)))
legend(x="right", legend=c("low","  |","  |","  |","high"), fill=c("white",topo.colors(4)))


# Ordination with these clusters
my_distances <- c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")
par(mfrow=c(2,4),mar=c(2,2,5,1), mgp=c(0.5,0,0))
## Default parameters, all distances
for (m in my_distances) {
  if (m == "aitchison") d <- vegdist(mymat+0.000000001, method=m)
  else d <- vegdist(mymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       xlab="PC 1", ylab="PC 2", asp=1, xaxt="n", yaxt="n")
  title(m, adj = 0.5, line = 1)
}
mtext("Default parameters", side = 3, line = -2, outer = TRUE, cex=1.5)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0))


# Distribution of distances
par(mfrow=c(2,4),mar=c(3,3,5,1), mgp=c(1.5,0.5,0))
for (m in my_distances) {
  if (m == "aitchison") {
    d <- vegdist(mymat+0.000000001, method=m)
  } else {
    d <- vegdist(mymat, method=m)
  }
  plot(density(d), lwd=3, col="black", main="", xlab="Value", ylab="Freq")
  title(m, adj=0.5, line=1)
}
mtext("Distance Distributions (default params)", side = 3, line = -2, outer = TRUE, cex=1.5)


# What happens after we apply LGD?



# Scratch : Beta distribution plots
# --- Set up x values
# len_x <- seq(0,1,by=0.025)    # x values for the beta distributions
# --- Plotting the beta distribution with various values
# plot(len_x, dbeta(len_x, 0.5,1), type="l", col="blue", lwd=3, ylab="Prob Dist", xlab="X", main="PDF of Beta", ylim=c(0,3.5))
# text(0.17, 1.7, label="beta(0.5,1)", col="blue")
# lines(len_x, dbeta(len_x, 0.5,0.5), col="red", lwd=3)
# text(0.58, 0.5, label="beta(0.5,0.5)", col="red")
# lines(len_x, dbeta(len_x, 1,1), col="purple", lwd=3)
# text(0.82, 1.2, label="beta(1,1)", col="purple")
# lines(len_x, dbeta(len_x, 3,3), col="darkgreen", lwd=3)
# text(0.61, 1.95, label="beta(3,3)", col="darkgreen")
# lines(len_x, dbeta(len_x, 3,1.5), col="cyan2", lwd=3)
# text(0.8, 2, label="beta(3,1.5)", col="cyan2")
# lines(len_x, dbeta(len_x, 6,8), col="orange", lwd=3)
# text(0.55, 2.75, label="beta(6,8)", col="orange")
# lines(len_x, dbeta(len_x, 3,10), lwd=3, lty="dashed")
# text(0.32, 3.2, label="beta(3,10)")
# --- Output of the the beta functions
# par(mar=c(4,2,2,2), mfrow=c(2,2))
# plot(dbeta(len_x, 6,8), main="Density: dbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(pbeta(len_x, 6,8), main="Distribution: pbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(qbeta(len_x, 6,8), main="Quantile Func: qbeta(6,8)", pch=16, cex=2, xlab="", ylab="")
# plot(density(rbeta(len_x, 6,8)), main="Random Deviates: rbeta(6,8)", lwd=3, xlab="", ylab="")
# par(mar=c(5.1,4.1,4.1,2.1), mfrow=c(1,1))

