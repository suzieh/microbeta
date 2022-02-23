# pca_simulation.r
# Principal Component Analysis on Simulation Data
# Knights Lab - University of Minnesota
# July 2019
# usage : pca_simulation.r

##### Set Up #####
#library(optparse)
library(lle)
library(ggplot2)
library(ggpubr)


##### Parse Command Line #####
# option_list <- list(make_option(c("-i", "--input_table"), type="character",
#                                 help="Path to input file. Expects otu table in tab-delimited txt file format.") )
# opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))
dat_all <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun.txt", sep="\t", header=T)
dat <- dat_all[1:200,] # remove noise for now
colors <- viridis(nrow(dat))
names(colors) <- rownames(dat)


##### Helpful Functions #####
# colorfunc <- colorRampPalette(c("white", "navy"))
# calc.perc.var <- function (eigen, dimension) {
#   percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
#   return(percents[dimension])
# }
plot.pc <- function (pc_cmd, grad_len) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  p <- ggplot(pc_cmd, aes(x=Dim1, y=Dim2, color=factor(color))) +
    geom_point(size = 3) +
    scale_color_manual(values = unname(colors[rownames(pc_cmd)])) +
    labs(title = paste0("ends at ", grad_len*5), x="", y="") +
    theme_classic() +
    theme(legend.position = 'none')
  return(p)
}
plot.pc1 <- function (pc_cmd) {
  pc_cmd <- as.data.frame(pc_cmd)
  colnames(pc_cmd)[1:2] <- c("Dim1", "Dim2")
  pc_cmd$color <- 1:nrow(pc_cmd)
  pc_cmd$sample <- as.numeric(gsub("sample.", "", rownames(pc_cmd)))
  p <- ggplot(pc_cmd, aes(x=sample, y=Dim1, color=factor(color))) +
    geom_point(size = 3) +
    scale_color_manual(values = unname(colors[rownames(pc_cmd)])) +
    labs(title = "", x = "", y = "") + 
    theme_classic() + 
    theme(legend.position = 'none')
  return(p)
}


##### Run Principal Components Analysis #####
# Compute Prinicpal Components with raw data
pca_out <- prcomp(dat, center=T, scale.=F) # centered at 0, variables scaled to have unit variance 
screeplot(pca_out) # scree plot (no elbow due to homogenity of data)
# Compute Principal Components in pieces of gradient
plot_list <- list() #change plot output dims: par(mfrow=c(4,4), oma=c(0,0,0,0), mar=c(2,2,2,2))
for(grad_len in seq(25, nrow(dat), by=25)) {
  n <- (grad_len/25) + ((grad_len-25)/25)
  inds <- seq.int(1, grad_len, length.out=25) # set indices to plot
  bc <- vegdist(dat[inds,]) # bray curtis distances
  pc_cmd <- cmdscale(bc)    # principal components analysis of bray curtis distances
  plot_list[[n]] <- plot.pc(pc_cmd, grad_len) # basic plot of PCs : plot(pc_cmd, main=paste0("ends at ", grad_len))
  #plot(as.numeric(dist(pc_cmd)),as.numeric(bc))   # plot PC distances vs. Bray Curtis distances 
  plot_list[[n+1]] <- plot.pc1(pc_cmd) # plot just PC1: plot(inds, pc_cmd[,1], main = "PC 1")
}
grid.arrange(grobs=plot_list, ncol=4, top=text_grob("Growing Gradient Length", face='bold')) #par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))


##### Extra Analyses - Bray Curtis and Jaccard #####
# plot(as.numeric(dist(pc_cmd)),as.numeric(bc), main="Comparing Actual vs. PC", xlab="actual distance", ylab="bray distance")
# # Look at Difference b/w Bray-Curtis & Jaccard (still omitting noise)
# pc_bc <- cmdscale(vegdist(dat))
# pc_jac <- cmdscale(vegdist(dat, method="jaccard")) # just curious about jaccard
# plot(pc_bc, main = "Bray-Curtis Dissimilarity", xlab = "PC1", ylab= "PC2")
# plot(pc_jac, main = "Jaccard Index of Bray-Curtis Dissimilarity", xlab = "PC1", ylab= "PC2")
# # nice version of Bray-Curtis for all samples (incorporating noise)
# pc_bc_all <- cmdscale(vegdist(dat_all))
# ggplot(data.frame(X=pc_bc_all[,1],Y=pc_bc_all[,2], color=rep(1:200,2), noise=c(rep("no",200), rep("yes",200))),
#        aes(x=X, y=Y, color=factor(color), shape=factor(noise), size=factor(noise))) +
#   geom_point(alpha=0.5) +
#   scale_color_viridis_d(nrow(pc_bc)/2) +
#   scale_shape_manual(values=c(1, 18)) +
#   scale_size_manual(values=c(15, 10)) +
#   labs(title="PCA : Bray Curtis Dissimilarity",
#        x="PC1",
#        y="PC2",
#        shape="noise") +
#   guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
#   theme_classic() + NULL
# # Could run this again on CLR transformed data, but seems like a waste (discussed in meeting in late August)
# 
