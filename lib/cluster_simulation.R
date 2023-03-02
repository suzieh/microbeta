# Generating a group-wise clusters simulation
# Suzie Hoops
# Last Updated: March 2023

library(stats)
library(plot.matrix)
library(vegan)
library(scales)
set.seed(125)

##### Function to create clusters dataset #####
# Create nclust clusters from nsamp samples
## note: balanced group sizes, but could add another parameter for inbalanced groups
create_clusters <- function (nsamp=90, nclust=3, min_n_otu=150, overlap=0, beta_ab=c(4,8), random_deviates=F) {
  # Initialize sample sizes
  if (overlap > 0 & overlap < 1) {
    notu_overlap <- min_n_otu * overlap                 # maximum number of overlapping OTUs
    notu_each <- max(c(50, round(min_n_otu / nclust)))  # approx. number of otus each sample
    notu_each <- notu_each + round((notu_overlap * (nclust-1)) / nclust) # adjust for overlap
    ncol_otu <- (notu_each * nclust) - (notu_overlap * (nclust-1))
  } else if (overlap > 1) {
    stop("overlap parameter must be a value between 0 and 1.")
  } else {
    notu_overlap <- 0                                   # no overlap case
    notu_each <- max(c(50, round(min_n_otu / nclust)))  # approx. number of otus each sample
    ncol_otu <- notu_each * nclust
  }
  # Initialize data matrix
  len_x <- seq(0, 1, length.out=1000 )                  # x values for the beta function
  matrix <- matrix(data=0, nrow=nsamp, ncol=ncol_otu)   # empty matrix to fill with samples
  sample_breaks <- c(rep(1:nclust, each=nsamp/nclust), rep(nclust, nsamp %% nclust))
  col_breaks <- lapply(seq(1,ncol_otu, by=notu_each - notu_overlap), function (start) {seq(start, start+notu_each-1, by=1)})
  cols <- lapply(sample_breaks, function (x) {col_breaks[[x]]})
  # Sampling: beta distribution to approximate a real microbiome dataset
  for (i in 1:nsamp) {
    if (random_deviates == T) {
      pdf <- density(rbeta(len_x, beta_ab[1], beta_ab[2]))$y  # create PDF from beta distribution
    } else {
      pdf <- dbeta(len_x, beta_ab[1], beta_ab[2])             # create PDF from beta distribution
    }
    curr_samp <- pdf[sort(runif(notu_each, 1, length(pdf)))]  # random sample from PDF
    curr_samp[curr_samp < 0.001] <- 0                         # get rid of low abundance OTUS
    curr_samp <- curr_samp / sum(curr_samp)                   # normalize sample
    matrix[ i, cols[[i]] ] <- curr_samp                       # add the sample to out matrix
  }
  # Clean up matrix and return
  matrix <- matrix[rowSums(matrix) > 0,]  # remove empty rows (should not happen, but safety check)
  matrix <- matrix[,colSums(matrix) > 0]  # remove empty columns
  return(matrix)
}
# Plotting output : no overlap
par(mar=c(5.1,4.1,4.1,4.1))
mymat <- create_clusters(random_deviates=F)
plot(mymat, col=c("yellow","green","cyan","royalblue","blue"), breaks=seq(0.000001, max(mymat), length.out=6),
     main="Clusters (n=3) : beta(4,8)", border=NA, axis.col=list(cex.axis=0.8), axis.row=list(cex.axis=0.8),
     key=list(cex.axis=0.8), fmt.key="%.3f")
rmymat <- create_clusters(random_deviates=T)
plot(rmymat, col=c("yellow","green","cyan","royalblue","blue"), breaks=seq(0.000001, max(rmymat), length.out=6),
     main="Clusters (n=3) : rbeta(4,8) - Random Dev.", border=NA, axis.col=list(cex.axis=0.8), axis.row=list(cex.axis=0.8),
     key=list(cex.axis=0.8), fmt.key="%.3f")
# Plotting output : ~ 10% overlap
tmymat <- create_clusters(overlap=0.1)
plot(tmymat, col=c("yellow","green","cyan","royalblue","blue"), breaks=seq(0.000001, max(tmymat), length.out=6),
     main="Clusters (n=3) : beta(4,8) - 10% Overlap", border=NA, axis.col=list(cex.axis=0.8), axis.row=list(cex.axis=0.8),
     key=list(cex.axis=0.8), fmt.key="%.3f")
# Plotting output : ~ 30% overlap
omymat <- create_clusters(overlap=0.3)
plot(omymat, col=c("yellow","green","cyan","royalblue","blue"), breaks=seq(0.000001, max(omymat), length.out=6),
     main="Clusters (n=3) : beta(4,8) - 30% Overlap", border=NA, axis.col=list(cex.axis=0.8), axis.row=list(cex.axis=0.8),
     key=list(cex.axis=0.8), fmt.key="%.3f")
# Plotting output : ~50% overlap
o2mymat <- create_clusters(overlap=0.5)
plot(o2mymat, col=c("yellow","green","cyan","royalblue","blue"), breaks=seq(0.000001, max(o2mymat), length.out=6),
     main="Clusters (n=3) : beta(4,8) - 50% Overlap", border=NA, axis.col=list(cex.axis=0.8), axis.row=list(cex.axis=0.8),
     key=list(cex.axis=0.8), fmt.key="%.3f")


# Ordination with these clusters
par(mfrow=c(2,4),mar=c(4,2,2,2), mgp=c(0.5,0,0))
## Defaults
for (m in c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")) {
  if (m == "aitchison") d <- vegdist(mymat+0.000000001, method=m)
  else d <- vegdist(mymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2.5, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       main=m, xlab="PC 1", ylab="PC 2", asp=3)
}
mtext("Default parameters", side = 3, line = -21, outer = TRUE)
## Random Deviates
for (m in c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")) {
  if (m == "aitchison") d <- vegdist(rmymat+0.000000001, method=m)
  else d <- vegdist(rmymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2.5, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       main=m, xlab="PC 1", ylab="PC 2", asp=T)
}
## Overlap 10%
for (m in c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")) {
  if (m == "aitchison") d <- vegdist(tmymat+0.000000001, method=m)
  else d <- vegdist(tmymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2.5, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       main=m, xlab="PC 1", ylab="PC 2", asp=T)
}
## Overlap 30%
for (m in c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")) {
  if (m == "aitchison") d <- vegdist(omymat+0.000000001, method=m)
  else d <- vegdist(omymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2.5, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       main=m, xlab="PC 1", ylab="PC 2", asp=T)
}
## Overlap 50%
for (m in c("bray", "jaccard", "aitchison", "euclidean", "manhattan", "gower", "chisq", "canberra")) {
  if (m == "aitchison") d <- vegdist(o2mymat+0.000000001, method=m)
  else d <- vegdist(o2mymat, method=m)
  pc <- cmdscale(d, k=2)
  plot(pc[,1], pc[,2], pch=16, cex=2.5, col=alpha(rep(c("purple", "orange", "blue"), each=30),0.6),
       main=m, xlab="PC 1", ylab="PC 2", asp=T)
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0))



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
