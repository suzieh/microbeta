
# Scratch: Reconstruction error for PCA
## note: impossible for PCoA... unless we try reconstructing the distance matrix.
## note: might be erroneous even for PCA, see "NONASYMPTOTIC UPPER BOUNDS FOR
##       THE RECONSTRUCTION ERROR OF PCA" by Reiss & Wahl
## compute a traditional PCA
pca = prcomp(t(soil_norm))
plot(pca$x[,1], pca$x[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs, main="PCA", xlab="PC 1", ylab="PC 2")
## compute a traditional SVD
soil_centered = sweep(t(soil_norm), 2, colMeans(soil_norm), "-")
sv <- svd(soil_centered)
plot(sv$u[,1], sv$u[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs, main="SVD", xlab="U 1", ylab="U 2")
## compute PCA how we would normally - matches the above!
pc = cmdscale(vegdist(t(soil_norm), "euclidean"), k=80, eig=T)
plot(pc$points[,1]*-1, pc$points[,2], pch=soil_shps, col=soil_cols, bg=soil_bgs, main="PCA my way", xlab="PC 1", ylab="PC 2")
## compute PCA reconstruction error:
projection <- ( pca$x %*% t(pca$rotation) ) - mean(as.matrix(soil_norm))
loss_matrix <- (t(soil_norm) - projection)^2
mse <- mean(loss_matrix)
mse 


# Scratch expand reconstruction error to PCoA???
## what if we compute a pcoa (not euclidean)
dist <- vegdist(t(soil_norm), "bray")
pcoa <- cmdscale(dist, k=80, eig=T)
## check that our symmetric covariance matrix is diagonalizable
v <- eigen(as.matrix(dist))$vectors
d <- diag(eigen(as.matrix(dist))$values)
out <- v %*% d %*% Matrix::solve(v) # dot product to get 88 x 88 covariance matrix back
out <- round(out, 7)
sum(abs(out - as.matrix(dist)) < 1e-6) # essentially identical matrices, so yes diagonizable
## can we create a projection from this v and d? or some variation?
## this v and d don't have the original dimensions... but can we reconstruct the distance matrix instead?
C = t(soil_norm) %*% as.matrix(soil_norm) # 88 x 88 covariance matrix
v <- eigen(C)$vectors
d <- diag(eigen(C)$values)
out <- v %*% d %*% Matrix::solve(v)
sum(abs(out - C) < 1e-6) # basically identical
x <- v %*% sqrt(d)