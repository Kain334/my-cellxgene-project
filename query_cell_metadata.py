import cellxgene_census

# Open the census database
with cellxgene_census.open_soma() as census:
    
    # Read the cell metadata as a slice, filtering female cells of specific types
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )
    
    # Concatenate results into a pyarrow Table
    cell_metadata = cell_metadata.concat()
    
    # Convert the result to a pandas DataFrame
    cell_metadata = cell_metadata.to_pandas()

    # Print the DataFrame
    print(cell_metadata)

