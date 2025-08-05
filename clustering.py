import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
import argparse

# === Argument parsing === #
parser = argparse.ArgumentParser(
    description="Cluster UMAP 3D projections using K-Means and plot the results."
)
parser.add_argument(
    "-k", "--clusters", type=int, required=True,
    help="Number of clusters to generate with K-Means (acceptable range: 2–100)"
)
args = parser.parse_args()

if not (2 <= args.clusters <= 100):
    raise ValueError("Number of clusters (-k) must be between 2 and 100.")

# === Load UMAP 3D projections from bio_embeddings === #
df = pd.read_csv("./umap_projections/projected_embeddings_file.csv")

df.rename(columns={
    "component_0": "x",
    "component_1": "y",
    "component_2": "z"
}, inplace=True)

# Sanity check
df.dropna(subset=["x", "y", "z"], inplace=True)
if df.empty:
    raise ValueError("No valid UMAP coordinates found. Check your projection file.")

# === K-Means Clustering === #
kmeans = KMeans(n_clusters=args.clusters, random_state=42)
df['cluster'] = kmeans.fit_predict(df[['x', 'y', 'z']])

# Save clustered output
df.to_csv("clustered_umap_3D_kmeans.csv", index=False)

# === Plotting === #
unique_clusters = np.sort(df['cluster'].unique())
num_clusters = len(unique_clusters)
cluster_map = {k: i for i, k in enumerate(unique_clusters)}
df['cluster_index'] = df['cluster'].map(cluster_map)

def get_spaced_colors(n, cmap_name='nipy_spectral', min_gap=3):
    cmap = plt.get_cmap(cmap_name)
    indices = np.linspace(0, 1, n * min_gap)[::min_gap]
    return cmap(indices[:n])

colors = get_spaced_colors(num_clusters)
cmap = ListedColormap(colors)
norm = BoundaryNorm(np.arange(-0.5, num_clusters + 0.5, 1), cmap.N)

fig = plt.figure(figsize=(12, 10), dpi=1200)
ax = fig.add_subplot(111, projection='3d')

scatter = ax.scatter(
    df['x'], df['y'], df['z'],
    c=df['cluster_index'],
    cmap=cmap,
    norm=norm,
    s=20,
    alpha=0.9
)

ax.view_init(elev=20, azim=135)
ax.set_title(f"3D UMAP + K-Means Clustering (k={args.clusters})", fontsize=14, fontweight='bold')
ax.set_xlabel("UMAP-1", fontsize=12)
ax.set_ylabel("UMAP-2", fontsize=12)
ax.set_zlabel("UMAP-3", fontsize=12)

ax.grid(False)
for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    axis.pane.fill = False

# Colorbar
cbar = plt.colorbar(scatter, ticks=np.arange(num_clusters), shrink=0.6, pad=0.1)
cbar.set_label('Cluster ID', fontsize=10)
cbar.ax.set_yticklabels(unique_clusters.astype(str))
cbar.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.savefig("clustered_umap_3D_kmeans.png", format='png', dpi=1200, bbox_inches='tight')

