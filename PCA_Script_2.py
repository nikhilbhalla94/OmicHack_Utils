import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Perform PCA analysis on gene expression data.')
parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file containing the data matrix.')
args = parser.parse_args()

# Load the data with low_memory=False
data = pd.read_csv(args.input, sep=',', header=None, low_memory=False)  # Explicitly set for CSV

# Extract sample IDs and group information
sample_ids = data.iloc[0, 1:].tolist()
groups = data.iloc[1, 1:].tolist()

# Extract gene expression data
gene_names = data.iloc[2:, 0].tolist()
expression_data = data.iloc[2:, 1:].astype(float)

# Perform PCA analysis
pca = PCA(n_components=14)  # Choose number of components
principal_components = pca.fit_transform(expression_data.T)  # Transpose to make samples rows

# Convert PCA results to a DataFrame
pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(14)], index=sample_ids)

# Save PCA coordinates to a CSV file
pca_df.to_csv('pca_coordinates.csv')

# Save explained variance to a CSV file
explained_variance = pd.DataFrame({'Principal Component': [f'PC{i+1}' for i in range(14)], 
                                   'Explained Variance': pca.explained_variance_ratio_})
explained_variance.to_csv('pca_variances.csv', index=False)

# Plotting the PCA
plt.figure(figsize=(7, 7))

# Sort groups for legend
sorted_groups = sorted(set(groups))

# Plot each group with adjusted data point size, border, and transparency
for group in sorted_groups:
    group_indices = [i for i, x in enumerate(groups) if x == group]
    plt.scatter(pca_df.iloc[group_indices, 0], pca_df.iloc[group_indices, 1], 
                label=group, 
                s=25,                 # Size of the data points
                edgecolor='black',    # Color of the border
                linewidth=0.5)       # Thickness of the border (adjust as needed)

# Calculate variance percentage for PC1 and PC2
pc1_variance = pca.explained_variance_ratio_[0] * 100
pc2_variance = pca.explained_variance_ratio_[1] * 100

# Set axis labels with variance percentages and use Arial font
plt.xlabel(f'PC1 ({pc1_variance:.2f}%)', fontsize=12, fontname='Arial')
plt.ylabel(f'PC2 ({pc2_variance:.2f}%)', fontsize=12, fontname='Arial')


# Move legend to the bottom left, sorted and remove the grid
plt.legend(loc='lower left', fontsize=10, frameon=False)
plt.grid(False)  # Remove grid lines

# Adjust border line thickness
ax = plt.gca()  # Get current axes
for spine in ax.spines.values():
    spine.set_linewidth(2)  # Set the thickness of the border lines (adjust as needed)

# Adjust ticks thickness
ax.tick_params(axis='both', which='both', width=2)  # Set the thickness of the ticks (adjust as needed)

plt.savefig('pca_plot.png')
plt.show()