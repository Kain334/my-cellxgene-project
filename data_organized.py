import urllib.request
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Download the dataset
urllib.request.urlretrieve('https://datasets.cellxgene.cziscience.com/5d871206-9489-4d9f-8106-94305ccb1c3a.h5ad', 'dataset.h5ad')

# Load the dataset
adata = sc.read_h5ad('dataset.h5ad')
print(adata)

# Check the list of obs and var attributes
obs_attributes = adata.obs.columns.tolist()
var_attributes = adata.var.columns.tolist()

# Slice female data
adata_slice_female = adata[adata.obs['sex'] == 'female', :]

# Set age stages for analysis
age_stages = [
    '57-year-old stage', '61-year-old stage', '62-year-old stage', '66-year-old stage', 
    '67-year-old stage', '71-year-old stage', '73-year-old stage', '75-year-old stage', 
    '77-year-old stage', '79-year-old stage', '81-year-old stage', '87-year-old stage', 
    '89-year-old stage'
]
filtered_adata = adata[adata.obs['development_stage'].isin(age_stages)].copy()
print("Filtered data count:", filtered_adata.shape)

# Divide data into younger and older age groups
young_group = filtered_adata[filtered_adata.obs['development_stage'].isin(
    ['57-year-old stage', '61-year-old stage', '62-year-old stage', '66-year-old stage', '67-year-old stage', '70-year-old stage']
)]
older_group = filtered_adata[filtered_adata.obs['development_stage'].isin(
    ['71-year-old stage', '73-year-old stage', '75-year-old stage', '77-year-old stage', '79-year-old stage', '81-year-old stage', '87-year-old stage', '89-year-old stage']
)]

# Define attributes to compare and compute comparison
attributes_to_compare = ['nCount_RNA', 'nFeature_RNA', 'percent.mt']
comparison = pd.DataFrame({
    'Attribute': attributes_to_compare,
    'Young Group Mean': [young_group.obs[attr].mean() for attr in attributes_to_compare],
    'Older Group Mean': [older_group.obs[attr].mean() for attr in attributes_to_compare]
})

# Display comparison results
print("Comparison of Attributes Between Younger and Older Groups:")
print(comparison)
