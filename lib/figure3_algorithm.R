# lgd_source.r
# Figure : Algorithm visually explained
# last updated March 2023

library(igraph)
library(plotrix)
library(vegan)
source('lib/lgd_source.r')
set.seed(125)

# Small set of simulated data points
mini_sim <-read.table("data/old_gradient/fake_rel_abun_long_n100.txt", row=1, header=T, sep="\t")
mini_sim <- mini_sim[c(1,12,19,31,43,50,61,72,80,92,100),]
mini_cols <- viridis::viridis(11, alpha=0.8)


# Before & after plots
mini_d <- vegdist(mini_sim, method="euclidean")
## Before plot
mini_pc <- cmdscale(mini_d, k=2, eig=F)
#tiff("figures/figure3_algorithm/1_before_pcoa.tiff", width=2, height=2, units="in", res=800) # needs to be 1200 dpi for final
par(mgp=c(1,0,0), mar=c(2,2,1,1))
plot(mini_pc[,1],mini_pc[,2], xaxt='n', yaxt='n',
     xlim=range(mini_pc)+c(-0.02,0.02), ylim=range(mini_pc)+c(-0.01,0.02),
     col=mini_cols, pch=16, cex=6, asp=1, 
     xlab="PC 1", ylab="PC 2", cex.lab=2)
#dev.off()
## After plot
mini_lgd <- lg.dist(mini_d, neighborhood.radius=0.18)
mini_pc_lgd <- cmdscale(mini_lgd, k=2, eig=F)
plot(mini_pc_lgd[,1],mini_pc_lgd[,2], xaxt='n', yaxt='n',
     xlim=range(mini_pc_lgd)+c(-0.02,0.02), ylim=range(mini_pc_lgd)+c(-0.02,0.02),
     col=mini_cols, pch=16, cex=6, asp=1,
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)


# iGraph plots of distance space and the connections made
## Make graph object from the points
# NOTE: using different radius then best radius for a more informative example
adj <- matrix(0, nrow(mini_sim), nrow(mini_sim))
for(i in 1:nrow(adj)) {
  adj[i,which(as.matrix(mini_d)[i,] <= 0.24823)] <- 1 
}
adj[adj>0] <- as.matrix(mini_d)[adj>0]
g <- graph.adjacency(adj, weighted=T, mode='undirected')
E(g)

## Prepare edges for one node
edge_colors3 <- rep("white",22)
edge_colors3[c(7,8,10:12)] <- "darkgrey"
edge_widths3 <- rep(0,22); edge_widths3[c(7,8,10:12)] <-3

## Plot graph and add radius circle if applicable
### note: change the edge.color parameter and re-run full script
###       to make sure point placement is the same. Run with "white"(for part2),
###       "darkgrey"(for part4) and edge_colors3 (for part3), edge_colors5 (for part5).
###       edge.width should be 3 or one of the edge_width variables.
plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
     edge.color=edge_colors3, edge.width=edge_widths3)
draw.circle(x=0.02, y=-0.07, radius=0.7, lwd=3, lty="dashed")
text(x=0.5, y=-0.82, labels="radius = 0.25", cex=0.9)

## Plot paths between node 2 and 11 (commented out if not applicable)
# edge_colors5 <- rep("white",22)
# edge_colors5[c(4,7,12,19)] <- "darkgrey"
# edge_widths5 <- rep(0,22); edge_widths5[c(4,7,12,19)] <- 3
# print(paste0("original dist: ",round(as.matrix(mini_d)[2,11],3))) # original: 0.284
# print(paste0("new dist: ",round(shortest.paths(g)[2,11],3)))      # LGD: 0.783
# plot(g, vertex.label=NA, vertex.color=mini_cols, vertex.frame.color=mini_cols, vertex.size=22,
#      edge.color=edge_colors5, edge.width=edge_widths5)
# lines(x=c(-0.9,0.71), y=c(-0.75,0.9), lwd=3, lty="dashed", col="black")
# legend("bottomright", legend=c("original = 0.284", "LGD = 0.783"), title="Inter-node Distance",
#        col=c("black","darkgrey"), lty=c("dashed","solid"), lwd=2, bty="n", cex=0.9)





