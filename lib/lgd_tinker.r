# lgd_tinker.r
# Tinkering with Local Gradient Distance
# Knights Lab - University of Minnesota
# September 2019
# usage : source('lgd_tinker.r')

##### Hanging Questions #####
# What's with the weights?
# Fix neighborhood.size vs. neighborhood radius
# Add Error message for NULL result
# Add warning for number of links created
# Add warning for links added > twice length of local edge
# calculate based on normalized table or raw counts? What's the difference?


##### Set Up #####
library(vegan)
library(ggplot2)
require('igraph', warn.conflicts=FALSE, quietly=TRUE) # load package quietly
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun_long.txt", sep="\t", header=T)
dat <- round(dat[c(5,10,15,20, 45,50,52,55, 75,82,88,90,92,96),]*10000) # select interesting set, convert to counts
rownames(dat) <- paste0("sample.", 1:nrow(dat))


##### N-based LGD #####
# Parameters
d <- vegdist(dat)
neighborhood.size = 3
neighborhood.radius = NULL
weighted = TRUE

# Function Contents
# step 1: neighborhood settings : uses radius if both provided, ranges allowed
if(is.null(neighborhood.size) && is.null(neighborhood.radius)){
  neighborhood.sizes <- 3:10 # default try n=3:10
} else if (is.null(neighborhood.radius)) {
  neighborhood.sizes <- neighborhood.size
  use.r <- FALSE
} else {
  neighborhood.sizes <- neighborhood.radius
  use.r <- TRUE
}
# WITH CURRENT PARAMS: n.ss = 3, use.r = F

# step 2/3: create graph, check validity (components)
#   calculate number of eigenvalues == 0 in graph;
#   this is the number of connected components.
ix <- 1
is.valid <- FALSE
while(ix <= length(neighborhood.sizes) && !is.valid){
  ns <- neighborhood.sizes[ix]
  # step 2: create graph g <- lg.graph(d, neighborhood.size=ns, use.r= use.r, weighted=weighted)
  neighborhood.size <- ns
  use.r <- use.r
  weighted <- weighted
  # start of lg.graph
  d <- as.matrix(d)
  adj <- matrix(0, nrow(d), nrow(d))
  for(i in 1:nrow(adj)) {
    if(!use.r) { adj[i,order(d[i,])[1:(neighborhood.size+1)]] <- 1 }
    else adj[i,which(d[i,] < neighborhood.size)] <- 1
  }
  if(weighted) adj[adj>0] <- d[adj>0] else weighted <- NULL
  rownames(adj) <- rownames(d); colnames(adj) <- colnames(d);
  g <- graph.adjacency(adj, weighted=weighted, mode='undirected')
  # end of lg.graph
  eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
  is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1
  ix <- ix + 1
}

# step 4: greedily add bridge if invalid. Need to throw a
# warning if bridge added is more than twice the length
# of a local edge (!)
if (!is.valid) {
  # start of greedy.connect
  dc <- split(1:nrow(as.matrix(d)), components(g)$membership)
  mstree <- graph.adjacency(mst(d), weighted=TRUE, mode='undirected')
  me <- split(as_edgelist(mstree), seq(gsize(mstree))); ge <- split(as_edgelist(g), seq(gsize(g)));
  new_edges <- unname(unlist(me[!(me %in% ge)]))
  w <- d[sapply(seq(1,length(new_edges),2), function (i) return((which(rownames(d) == new_edges[i])-1)*nrow(d)+which(colnames(d) == new_edges[i+1])))]
  warning(paste0("Disconnected graph with given neighborhood size. Adding ", length(new_edges)/2, " edges."), call.=FALSE)
  if (any(w > 2*max(E(g)$weight))) warning("Some added edge weights exceed twice length of neighborhood edges.")
  g <- add_edges(g, new_edges, attr = list("weight"=w))
  # end of greedy.connect
  eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
  is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1
}
# return distance matrix based on this graph
lgd <- shortest.paths(g)
# look at output
cat('Calculating local gradient distance...\n')
lgd <- lgd
cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd <- cmdscale(lgd)
cat('Plotting PCoA of original distances...\n')
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances", xlab="PC1", ylab="PC2")
cat('Plotting PCoA of transformed distances...\n')
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances", xlab="PC1", ylab="PC2")




##### Testing Function formats #####
"lg.dist" <- function(d, neighborhood.size=NULL, neighborhood.radius=NULL, weighted=TRUE) {
  # note: master function, calculates distances based on connected graph of neighborhoods
  # note: prefers to use radius, but if nothing provided tries n = 3:10
  if(is.null(neighborhood.size) && is.null(neighborhood.radius)){
    neighborhood.sizes <- 3:10 # default try n=3:10
  } else if (is.null(neighborhood.radius)) {
    neighborhood.sizes <- neighborhood.size
    use.r <- FALSE
  } else {
    neighborhood.sizes <- neighborhood.radius
    use.r <- TRUE
  }
  # note: try to create connected graph with provided neighborhood ranges.
  ix <- 1
  is.valid <- FALSE
  while(ix <= length(neighborhood.sizes) && !is.valid){
    ns <- neighborhood.sizes[ix]
    g <- lg.graph(d, neighborhood.size=ns, use.r= use.r, weighted=weighted)
    eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
    is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1
    ix <- ix + 1
  }
  # note: connect graph if still invalid
  if (!is.valid) {
    g <- greedy.connect(d, g, neighborhood.sizes[ix-1])
    eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
    is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1
  }
  # note: return dist object
  lgd <- shortest.paths(g)
  if (is.null(lgd)) {warning("Resulted in null graph."); return(NULL);}
  return(lgd)
}

"lg.graph" <- function(d, neighborhood.size=NULL, use.r=TRUE, weighted=TRUE) {
  # note: creates graph with provided neighborhood metric
  require('igraph', warn.conflicts=FALSE, quietly=TRUE)
  d <- as.matrix(d)
  adj <- matrix(0, nrow(d), nrow(d))
  for(i in 1:nrow(adj)) {
    if(!use.r) { adj[i,order(d[i,])[1:(neighborhood.size+1)]] <- 1 }
    else adj[i,which(d[i,] <= neighborhood.size)] <- 1
  }
  if(weighted) adj[adj>0] <- d[adj>0] else weighted <- NULL
  rownames(adj) <- rownames(d); colnames(adj) <- colnames(d);
  g <- graph.adjacency(adj, weighted=weighted, mode='undirected')
  return(g)
}

"greedy.connect" <- function (d, g, ns) {
  # note: greedily connects disconnected components with minimum distance between components
  # note: prints warning to output to notify user of this operation
  dc <- split(1:nrow(as.matrix(d)), components(g)$membership)
  mstree <- graph.adjacency(mst(d), weighted=TRUE, mode='undirected')
  me <- split(as_edgelist(mstree), seq(gsize(mstree))); ge <- split(as_edgelist(g), seq(gsize(g)));
  new_edges <- unname(unlist(me[!(me %in% ge)]))
  w <- d[sapply(seq(1,length(new_edges),2), function (i) return((which(rownames(d) == new_edges[i])-1)*nrow(d)+which(colnames(d) == new_edges[i+1])))]
  warning(paste0("Disconnected graph with given neighborhood size (", ns, "). Adding ", length(new_edges)/2, " edges."), call.=FALSE)
  if (any(w > 2*max(E(g)$weight))) warning("Some added edge weights exceed twice length of neighborhood edges.")
  g <- add_edges(g, new_edges, attr = list("weight"=w))
  return(g)
}

# TEST 1 : Disconnected, n=3
cat('Calculating local gradient distance...\n')
lgd <- lg.dist(d, neighborhood.size=3)
cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd <- cmdscale(lgd)
cat('Plotting PCoA of original distances...\n')
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances", xlab="PC1", ylab="PC2")
cat('Plotting PCoA of transformed distances...\n')
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances", xlab="PC1", ylab="PC2")

# TEST 2 : Disconnected, r=0.2:0.4
plot(density(d))
cat('Calculating local gradient distance...\n')
lgd <- lg.dist(d, neighborhood.radius=c(0.2, 0.3, 0.4))
cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd <- cmdscale(lgd)
cat('Plotting PCoA of original distances...\n')
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances", xlab="PC1", ylab="PC2")
cat('Plotting PCoA of transformed distances...\n')
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances", xlab="PC1", ylab="PC2")

# TEST 3 : Disconnected, r=minimium distance for connected groups (cheating based on known distances)
plot(density(d))
r <- max(d["sample.4", "sample.5"], d["sample.8", "sample.9"]) # get minimum r for connected graph
cat('Calculating local gradient distance...\n')
lgd <- lg.dist(d, neighborhood.radius=r)
cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')
pc.lgd <- cmdscale(lgd)
cat('Plotting PCoA of original distances...\n')
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances", xlab="PC1", ylab="PC2")
cat('Plotting PCoA of transformed distances...\n')
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances", xlab="PC1", ylab="PC2")




