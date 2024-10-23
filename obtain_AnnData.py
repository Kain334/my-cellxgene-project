import cellxgene_census

try:
    with cellxgene_census.open_soma(census_version="2024-07-01") as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
            obs_value_filter="sex == 'female' and cell_type in ['microglial cell', 'neuron']",
            obs_column_names=["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
        )
        print(adata)
except Exception as e:
    print(f"An error occurred: {e}")

