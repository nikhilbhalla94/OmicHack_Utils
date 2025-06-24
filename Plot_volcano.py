import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# Argument parser setup
parser = argparse.ArgumentParser(description="Generate a volcano plot with draggable labels.")
parser.add_argument("-i", "--input", required=True, help="Path to input CSV file")
parser.add_argument("-o", "--output", help="Output file name (SVG)")
parser.add_argument("-w", "--width", type=int, default=8, help="Width of the plot in inches")
parser.add_argument("-H", "--height", type=int, default=6, help="Height of the plot in inches")
parser.add_argument("--border", type=float, default=2.5, help="Border thickness (default: 2.5)")
parser.add_argument("--tickwidth", type=float, default=2, help="Tick width (default: 2)")
parser.add_argument("--ticklength", type=float, default=6, help="Tick length (default: 6)")
args = parser.parse_args()

# Read input file
data = pd.read_csv(args.input, index_col=0)

# Compute -log10(padj)
data["-log10(P_adj)"] = -np.log10(data["padj"])

# Drop rows with NaN values
data = data.dropna(subset=["log2FoldChange", "-log10(P_adj)"])

# Significance thresholds
pval_cutoff = 0.05
logfc_cutoff = 1.0

# Categorize genes
data["Significance"] = "Not Significant"
data.loc[(data["padj"] < pval_cutoff) & (data["log2FoldChange"] > logfc_cutoff), "Significance"] = "Upregulated"
data.loc[(data["padj"] < pval_cutoff) & (data["log2FoldChange"] < -logfc_cutoff), "Significance"] = "Downregulated"

# Define colors
color_map = {
    "Upregulated": "darkred",
    "Downregulated": "darkblue",
    "Not Significant": "#95A5A6"
}

# Output filename
output_file = args.output if args.output else os.path.splitext(args.input)[0] + ".svg"

# Create Matplotlib Figure
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(args.width, args.height))

# Scatter plot
ax.scatter(data["log2FoldChange"], data["-log10(P_adj)"], 
           c=data["Significance"].map(color_map), alpha=0.6)

# Add threshold lines
ax.axhline(-np.log10(pval_cutoff), color='grey', linestyle='dashed', lw=1.5, alpha=0.2)
ax.axvline(-logfc_cutoff, color='grey', linestyle='dashed', lw=1.5, alpha=0.2)
ax.axvline(logfc_cutoff, color='grey', linestyle='dashed', lw=1.5, alpha=0.2)

ax.set_xlabel(r"log$_{2}$ Fold Change", fontsize=14)
ax.set_ylabel(r"-log$_{10}$ P$_{adjusted}$", fontsize=14)


# Apply user-defined border thickness
for spine in ax.spines.values():
    spine.set_linewidth(args.border)  # Set border thickness

# Apply user-defined tick width and length
ax.tick_params(axis='both', which='both', width=args.tickwidth, length=args.ticklength)

# Store labels and connectors
texts = []
connectors = []

for idx, row in data.iterrows():
    if row["Significance"] in ["Upregulated", "Downregulated"]:
        ha = 'right' if row["log2FoldChange"] < 0 else 'left'
        color = 'darkred' if row["Significance"] == "Upregulated" else 'darkblue'

        # Create text label
        text = ax.text(row["log2FoldChange"], row["-log10(P_adj)"] + 0.5, 
                       idx, fontsize=8, ha=ha, va='center', color=color)

        # Create connector line
        line, = ax.plot([row["log2FoldChange"], row["log2FoldChange"]], 
                        [row["-log10(P_adj)"], row["-log10(P_adj)"] + 0.5], 
                        color=color, linestyle="dashed", lw=0.8, alpha=0.6)

        texts.append(text)
        connectors.append(line)

### **Enable Dragging for Labels & Update Connectors**
class DraggableAnnotations:
    def __init__(self, texts, connectors, data_points):
        self.texts = texts
        self.connectors = connectors
        self.data_points = data_points  # Store original data locations
        self.pressed = None
        self.cidpress = fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        for i, text in enumerate(self.texts):
            contains, _ = text.contains(event)
            if contains:
                self.pressed = i
                break

    def on_release(self, event):
        self.pressed = None
        fig.canvas.draw()

    def on_motion(self, event):
        if self.pressed is not None and event.xdata is not None and event.ydata is not None:
            text = self.texts[self.pressed]
            line = self.connectors[self.pressed]
            orig_x, orig_y = self.data_points[self.pressed]  # Keep one end fixed

            # Move label to new position
            text.set_position((event.xdata, event.ydata))

            # Update connector (fixed at data point)
            line.set_xdata([orig_x, event.xdata])
            line.set_ydata([orig_y, event.ydata])

            fig.canvas.draw()

# Activate draggable annotations
data_points = list(zip(data["log2FoldChange"], data["-log10(P_adj)"]))  # Store original data positions
dr = DraggableAnnotations(texts, connectors, data_points)

# Show interactive Matplotlib plot
plt.show()

# Save after interaction
plt.savefig(output_file, format="svg", dpi=1200, bbox_inches="tight")
print(f"Plot saved as: {output_file}")
