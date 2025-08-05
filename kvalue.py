import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="Input projected_embeddings_file.csv")
parser.add_argument("--kmin", type=int, default=2)
parser.add_argument("--kmax", type=int, default=40)
args = parser.parse_args()

df = pd.read_csv(args.input)
X = df[["component_0", "component_1", "component_2"]].values

inertia = []
k_range = range(args.kmin, args.kmax + 1)
for k in k_range:
    km = KMeans(n_clusters=k, n_init="auto", random_state=42)
    km.fit(X)
    inertia.append(km.inertia_)

# Save CSV
csv_out = "kmeans_elbow_values.csv"
pd.DataFrame({"k": list(k_range), "inertia": inertia}).to_csv(csv_out, index=False)

# Plot
plt.figure()
plt.plot(k_range, inertia, marker="o")
plt.xlabel("Number of clusters (k)")
plt.ylabel("Inertia")
plt.title("Elbow Method for Optimal k")
plt.grid(True)
plt.tight_layout()
plt.savefig("kmeans_elbow_plot.png", dpi=300)

