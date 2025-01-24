# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 12:26:18 2025

@author: lafields2
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rcParams['font.family'] = 'Times New Roman'
# Data (example)
data = {
    'Before De-N-Glycosylation': [60, 35, 30, 25, 59],
    'After De-N-Glycosylation': [98, 57, 39, 38, 91]
}
df = pd.DataFrame(data, index=['1', '2', '3', '4', '≥5'])

# Create a heatmap
plt.figure(figsize=(6, 4))
ax = sns.heatmap(
    df, 
    annot=True, 
    fmt='d', 
    cmap='flare',  # Change the color scheme here
    cbar=True, 
    linewidths=0.5,  # Add grid lines for better aesthetics
    linecolor='white'
)

# Move the x-axis labels to the top
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

# Rotate the y-axis labels vertically
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')

# Add title and labels
plt.title('Heatmap of De-N-Glycosylation', pad=35)
plt.xlabel('')
plt.ylabel('')

# Display the heatmap
plt.tight_layout()
plt.show()