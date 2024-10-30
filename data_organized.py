import urllib.request
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

# Plotting RNA count comparison
plt.figure(figsize=(8, 6))
sns.boxplot(data=adata.obs, x='disease', y='nCount_RNA')
plt.title("RNA Count Comparison Between Alzheimer and Normal Subjects")
plt.xlabel("disease")
plt.ylabel("nCount_RNA")
plt.show()

# Check the unique values in the 'sex' column
print(adata.obs['sex'].unique())

# Filter the dataset to focus on subjects categorized by 'sex'
male_group = adata[adata.obs['sex'] == 'male']
female_group = adata[adata.obs['sex'] == 'female']
male_counts = male_group.X.mean(axis=0)
female_counts = female_group.X.mean(axis=0)

# Create a DataFrame to store gene expression for easier comparison
genes = adata.var.index
gene_expression_comparison = pd.DataFrame({
    'Gene': genes,
    'Male': male_counts.A1 if hasattr(male_counts, "A1") else male_counts,
    'Female': female_counts.A1 if hasattr(female_counts, "A1") else female_counts
})

# Calculate the difference between male and female expression
gene_expression_comparison['Difference'] = gene_expression_comparison['Male'] - gene_expression_comparison['Female']
gene_expression_comparison_sorted = gene_expression_comparison.sort_values(by='Difference', ascending=False)

# Plot the top 10 genes with the largest difference in expression
top_genes = gene_expression_comparison_sorted.head(10)
plt.figure(figsize=(10, 6))
sns.barplot(x='Gene', y='Difference', data=top_genes)
plt.title('Top 10 Genes with Largest Expression Difference Between Male and Female Subjects')
plt.xticks(rotation=90)
plt.ylabel('Difference in Expression (Male - Female)')
plt.tight_layout()
plt.show()

# Distribution of nCount_RNA by Age Group
plt.figure(figsize=(8, 6))
plt.hist(young_group.obs['nCount_RNA'], bins=50, alpha=0.5, label='Young Group (<=70)')
plt.hist(older_group.obs['nCount_RNA'], bins=50, alpha=0.5, label='Older Group (>70)')
plt.title("Distribution of nCount_RNA by Age Group")
plt.xlabel("nCount_RNA")
plt.ylabel("Frequency")
plt.legend()
plt.show()

# Separate data by gender
adata_female = adata[adata.obs['sex'] == 'female'].copy()
adata_male = adata[adata.obs['sex'] == 'male'].copy()

# Calculate neighbors and UMAP
sc.pp.neighbors(adata_female)
sc.pp.neighbors(adata_male)
sc.tl.umap(adata_female)
sc.tl.umap(adata_male)

# Visualization - Female (Purple)
sc.pl.umap(adata_female, color='sex', title="UMAP - Female", color_map='Purples')

# Visualization - Male (Green)
sc.pl.umap(adata_male, color='sex', title="UMAP - Male", color_map='Greens')

