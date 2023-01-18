# lgd_source.r
# Local Gradient Distance
# Knights Lab - University of Minnesota
# Created August 2019
# usage : source('lgd_source.r')


##### Functions #####
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
    if (clusters(g)$no == 1) {
      is.valid <- TRUE # validity test
      # note: retired test for one eigenvalue = 0: sum(eigs < 10 * .Machine$double.eps) == 1
    }
    ix <- ix + 1
  }
  # note: connect graph if still invalid
  if (!is.valid) {
    g <- greedy.connect(d, g, neighborhood.sizes[ix-1])
    eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
    if (clusters(g)$no == 1) {
      is.valid <- TRUE # validity test
      # note: retired test for one eigenvalue = 0: sum(eigs < 10 * .Machine$double.eps) == 1
    }
    # previously
    # is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1   # validity test (one eigen val = 0)
    # if (!is.valid) {
    #   is.valid <- clusters(g)$no == 1                       # validity test (only one cluster, thus fully connected)
    # }
  }
  # note: return dist object
  lgd <- shortest.paths(g)
  # png(paste0("results/igraphs_gif/r_",round(neighborhood.radius,4),".png"), width=600, height=600)
  # plot(g, vertex.label=NA, vertex.color="purple", vertex.size=8,
  #      main=paste0("r = ", round(neighborhood.radius,4)))
  # dev.off()
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
  #dc <- split(1:nrow(as.matrix(d)), components(g)$membership)
  og <- graph.adjacency(as.matrix(d), weighted=T, mode='undirected')
  mstree <- mst(og) #graph.adjacency(mst(d), weighted=TRUE, mode='undirected')
  me <- split(as_edgelist(mstree), seq(gsize(mstree))); ge <- split(as_edgelist(g), seq(gsize(g)));
  new_edges <- unname(unlist(me[!(me %in% ge)]))
  d_idx <- sapply(seq(1,length(new_edges),2), function (i) {
    tmpd <- as.matrix(d)
    a <- which(rownames(tmpd) == new_edges[i])    # get row idx
    b <- which(colnames(tmpd) == new_edges[i+1])  # get column index
    return(c(a, b))
  })
  w <- sapply(1:ncol(d_idx), function (i) {as.matrix(d)[d_idx[1,i], d_idx[2,i]]})
  warning(paste0("Disconnected graph with given neighborhood size (", ns, "). Adding ", length(new_edges)/2, " edges."), call.=FALSE)
  if (any(w > 2*max(E(g)$weight))) warning("Some added edge weights exceed twice length of neighborhood edges.")
  g <- add_edges(g, new_edges, attr = list("weight"=w))
  return(g)
}

##### Usage #####
# cat('Calculating local gradient distance...\n')
# lgd <- lg.dist(distobj, neighborhood.size = 12, weighted=TRUE)
# cat('Calculating PCoA of original distances...\n')
# pc.d <- cmdscale(distobj)
# cat('Calculating PCoA of transformed distances...\n')
# pc.lgd <- cmdscale(lgd)
# cat('Plotting PCoA of original distances...\n')
# plot(pc.d, xlim=range(pc.d), ylim=range(pc.d), main="Original Distances")
# cat('Plotting PCoA of transformed distances...\n')
# plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd), main="Transformed Distances")

