# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 12:54:15 2025

@author: lafields2
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rcParams['font.family'] = 'Times New Roman'

lung_data = {
    'Before': [60, 35, 30, 25, 59],
    'After': [98, 57, 39, 38, 91]
}
heart_data = {
    'Before': [60, 35, 30, 25, 59],
    'After': [98, 57, 39, 38, 91]
}
kidney_data = {
    'Before': [60, 35, 30, 25, 59],
    'After': [98, 57, 39, 38, 91]
}
spleen_data = {
    'Before': [60, 35, 30, 25, 59],
    'After': [98, 57, 39, 38, 91]
}

# Create dataframes
lung_df = pd.DataFrame(lung_data, index=['1', '2', '3', '4', '≥5'])
heart_df = pd.DataFrame(heart_data, index=['1', '2', '3', '4', '≥5'])
kidney_df = pd.DataFrame(kidney_data, index=['1', '2', '3', '4', '≥5'])
spleen_df = pd.DataFrame(spleen_data, index=['1', '2', '3', '4', '≥5'])

# Create subplots
fig, axes = plt.subplots(1, 4, figsize=(7, 4), sharey=True)

# Plot each heatmap
sns.heatmap(lung_df, annot=True, fmt='d', cmap='flare', cbar=False, ax=axes[0])
axes[0].set_title('Lung')
axes[0].xaxis.set_ticks_position('top')
axes[0].xaxis.set_label_position('top')
axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0, ha='right')
axes[0].set_ylabel('PTM', fontsize=12)

sns.heatmap(heart_df, annot=True, fmt='d', cmap='flare', cbar=False, ax=axes[1])
axes[1].set_title('Heart')
axes[1].xaxis.set_ticks_position('top')
axes[1].xaxis.set_label_position('top')

sns.heatmap(kidney_df, annot=True, fmt='d', cmap='flare', cbar=False, ax=axes[2])
axes[2].set_title('Kidney')
axes[2].xaxis.set_ticks_position('top')
axes[2].xaxis.set_label_position('top')

sns.heatmap(spleen_df, annot=True, fmt='d', cmap='flare', cbar=True, ax=axes[3])
axes[3].set_title('Spleen')
axes[3].xaxis.set_ticks_position('top')
axes[3].xaxis.set_label_position('top')

fig.suptitle('Heatmap of De-N-Glycosylation', fontsize=16, y=1)

# Adjust layout
plt.tight_layout()
plt.show()
