import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Given ranges
min_value = 10
max_value = 40

# Create a figure and axis just for the color bar
fig, ax = plt.subplots(figsize=(0.1, 3))

# Set up the normalization for the color map
norm = mcolors.Normalize(vmin=min_value, vmax=max_value)

# Create a color bar without any associated plot
cbar = plt.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap="viridis"), cax=ax, orientation="vertical"
)

# Set the label for the color bar
cbar.set_label("Density", fontsize=15)

# Save the color bar to a PDF file
plt.savefig("colorbar_only.pdf", bbox_inches="tight")
