# Quick visualization of simulated dataset
idx <- rep(1:50, each=2)
idx <- idx + c(0,50)
heatmap(as.matrix(sim_norm[idx,]), Rowv=NA, Colv=NA, revC=F,
        col=c("white",viridis::viridis(49)))


