# lgd_source.r
# Local Gradient Distance
# Knights Lab - University of Minnesota
# Created August 2019, last updated March 2023
# usage : source('lgd_source.r')


##### Functions #####
"lg.dist" <- function(d, neighborhood.radius=NULL, weighted=TRUE, smooth=FALSE, epsilon=0.05) {
  # This master function calculates distances based on connected graph of neighborhoods
  # d : distance object to be adjusted
  # neighborhood.radius : radius to determine neighborhoods (NULL if algorithm should pick)
  # weighted : True/False use weighted edges in graph (recommended: TRUE)
  
  # Verify parameter validity
  ## distance object or distance matrix must be provided
  if (class(d) != "dist") {
    if (isSymmetric(as.matrix(d)) == TRUE) {
      d <- as.dist(d)
    } else {
      stop("'d' must be a distance object or symmetric distance matrix.")
    }
  }
  ## neighborhood radius must be a double or list of doubles or NULL (default)
  if (typeof(neighborhood.radius) != "double" & suppressWarnings(any(is.na(as.double(neighborhood.radius)) == TRUE))) {
    stop("'neighborhood.radius' must be numeric type 'double' or 'NULL' (default) to automatically find best neighborhood.")
  } else if (is.null(neighborhood.radius)) {
    neighborhood.sizes = seq(max(d), min(d), length.out=50) # try 50 radii in range of d
  } else {
    neighborhood.sizes = sort(as.double(neighborhood.radius), decreasing=T)
  }
  ## weighted should be True (default) or False
  if (typeof(weighted) != "logical") {
    stop("'weighted' must be logical type TRUE (default) or FALSE.")
  }

  #  Create connected graph with provided neighborhood range(s)
  ix <- 1         # index
  lgd <- NULL     # distances
  curr_corr <- -1 # optimization function score
  curr_ratio <- 1 # ratio degree:n
  curr_radius <- neighborhood.sizes[1]
  while(ix <= length(neighborhood.sizes)){
    ns <- neighborhood.sizes[ix]
    print(paste0("Trying radius = ", round(ns,3)))
    out <- lg.evaluate(d, neighborhood.size=ns, weighted=weighted) # list(lgd, is.valid, ratio, corr)
    print(paste0('   is.valid = ', out[[2]]))
    print(paste0('   ratio = ', out[[3]]))
    print(paste0('   corr = ', out[[4]]))
    if (!out[[2]]) {
      stop(paste0('Graph is invalid (r = ',ns,'), check all parameters are valid or try a different radius.'))
    }
    if ((out[[3]] >= 0.1) & (out[[4]] > (curr_corr + epsilon))) {
      print(paste0("   new best r found --> ", ns))
      lgd <- out[[1]]
      curr_ratio <- out[[3]]
      curr_corr <- out[[4]]
      curr_radius <- ns
    }
    if (out[[3]] < 0.1)
      break
    ix <- ix + 1
  }
  # Return dist object
  if (is.null(lgd)) {warning("Resulted in null graph."); return(NULL);}
  return(as.dist(lgd))
}

"lg.evaluate" <- function(d, neighborhood.size=NULL, weighted=TRUE) {
  # For a given radius value, returns graph & relevant inforamtion
  is.valid <- FALSE # default, must be proven otherwise
  n <- dim(as.matrix(d))[1] # number of samples
  g <- lg.graph(d, neighborhood.size=neighborhood.size, weighted=weighted)
  #eigs <- eigen(graph.laplacian(g), only.values=TRUE)$values
  if (clusters(g)$no == 1) {
    is.valid <- TRUE
  } else {
    # if invalid graph, use minimum spanning tree to add minimum number of edges
    g <- greedy.connect(d, g, neighborhood.size)
    if (clusters(g)$no == 1) {
      is.valid <- TRUE
    }
  }
  ratio <- mean(degree(g)) / n # record degree:n ratio
  if(is.valid) {
    lg_d <- shortest.paths(g)
    pc <- cmdscale(lg_d, k=3, eig=F)
    corr <- round(cor(dist(pc), as.dist(lg_d), method="pearson"),3) # use correlation b/w PC & LGD dists as optimization func
  } else {
    lg_d <- NULL
    corr <- -1
  }
  ## commented out : print graph
  # png(paste0("results/igraphs_gif/r_",round(neighborhood.radius,4),".png"), width=600, height=600)
  # plot(g, vertex.label=NA, vertex.color="purple", vertex.size=8,
  #      main=paste0("r = ", round(neighborhood.radius,4)))
  # dev.off()
  return(list(lg_d, is.valid, ratio, corr))
}

"lg.graph" <- function(d, neighborhood.size=NULL, weighted=TRUE) {
  # Creates graph with provided neighborhood metric
  require('igraph', warn.conflicts=FALSE, quietly=TRUE)
  d <- as.matrix(d)
  adj <- matrix(0, nrow(d), nrow(d))
  for(i in 1:nrow(adj)) {
    # retired neighborhood size n: adj[i,order(d[i,])[1:(neighborhood.size+1)]] <- 1
    adj[i,which(d[i,] <= neighborhood.size)] <- 1
  }
  if(weighted) adj[adj>0] <- d[adj>0] else weighted <- NULL
  rownames(adj) <- rownames(d); colnames(adj) <- colnames(d);
  g <- graph.adjacency(adj, weighted=weighted, mode='undirected')
  return(g)
}

"greedy.connect" <- function (d, g, ns) {
  # Greedily connects disconnected components with minimum distance between components
  # note: prints warning to output to notify user of this operation
  # retired method : dc <- split(1:nrow(as.matrix(d)), components(g)$membership)
  og <- graph.adjacency(as.matrix(d), weighted=T, mode='undirected')
  mstree <- mst(og) # retired: graph.adjacency(mst(d), weighted=TRUE, mode='undirected')
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

