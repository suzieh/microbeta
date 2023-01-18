import matplotlib.pyplot as plt
from sklearn import manifold, datasets
from matplotlib import cm
from matplotlib.colors import Normalize, to_hex
import numpy as np

# create dataset
sr_pts, sr_color = datasets.make_swiss_roll(300, random_state=0)

# visualize dataset
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
fig.add_axes(ax)
ax.scatter( sr_pts[:,0], sr_pts[:,1], sr_pts[:,2], c=sr_color, s=50, alpha=0.8 )
ax.set_title("Swiss Roll in Ambient Space")
ax.view_init(azim=-66, elev=12)
_ = ax.text2D(0.8, 0.05, s="n_samples=200", transform=ax.transAxes)
plt.show()

# export dataset (3D)
np.savetxt("Desktop/swiss_roll.txt", sr_pts)

# export colors
norm = Normalize(sr_color.min(), sr_color.max())
cmap = cm.viridis
sr_color_hex = []
for col in sr_color:
	sr_color_hex.append( to_hex( cmap(norm(col)) ) )
with open("Desktop/swiss_roll_colors.txt", "w") as color_file:
	color_file.write("\n".join(sr_color_hex))