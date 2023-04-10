# Playing around with Reconstruction Error in PCA
from sklearn.decomposition import PCA
import numpy as np
from numpy.testing import assert_array_almost_equal
import matplotlib.pyplot as plt

# Create a random dataset for this example
X_train = np.random.randn(100,50)
X_train.shape

# Fit PCA and compare to SVD (expected to be the same because PCA uses SVD here)
pca = PCA(n_components=30)  # initialize model and parameter of 30 dimensions
pca.fit(X_train)            # fit model to our dataset
U, S, VT = np.linalg.svd(X_train - X_train.mean(0), full_matrices=False)
assert_array_almost_equal(VT[:30], pca.components_)
## note: if above is not True, may need to shuffle U and V
# from sklearn.utils.extmath import svd_flip
# U, VT = svd_flip(U, VT)
# assert_array_almost_equal(VT[:30], pca.components_)

# Transform to the loadings in lower dimensional space
X_train_pca = pca.transform(X_train)  # typical PCA
X_train_pca.shape
X_train_pca2 = (X_train - pca.mean_).dot(pca.components_.T) # dot product with V approach
X_train_pca2.shape
assert_array_almost_equal(X_train_pca, X_train_pca2)
## visualize scree plot of explained variance
dims_n = np.arange(pca.n_components_) + 1
plt.plot(dims_n, pca.explained_variance_ratio_, 'o-', linewidth=2, color='blue')
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained')
plt.show()
## visualize first few dimensions
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(X_train_pca[:,1],X_train_pca[:,2],X_train_pca[:,3], marker='o')
ax.set_xlabel('PC 1')
ax.set_ylabel('PC 2')
ax.set_zlabel('PC 3')
plt.show()


# Project back to original dimensions (reconstruction)
X_projected = pca.inverse_transform(X_train_pca)
X_projected2 = X_train_pca.dot(pca.components_) + pca.mean_  # dot product for inverse
assert_array_almost_equal(X_projected, X_projected2)


# Calculate loss between original data and reconstruction (reconstruction error)
loss = np.sum((X_train - X_projected) ** 2, axis=1).mean() # mean squared error