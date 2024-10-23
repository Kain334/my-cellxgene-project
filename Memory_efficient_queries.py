import cellxgene_census
import tiledbsoma

# Open the census and query data for brain and male cells from human
with cellxgene_census.open_soma() as census:
    
    # Access human cell data
    human = census["census_data"]["homo_sapiens"]
    
    # Define a query for cells where tissue is brain and sex is male
    query = human.axis_query(
        measurement_name = "RNA",
        obs_query = tiledbsoma.AxisQuery(
            value_filter = "tissue == 'brain' and sex == 'male'"
        )
    )

    # Iterate over the query results
    iterator = query.X("raw").tables()
    
    # Example: process the first slice of the result
    raw_slice = next(iterator)
    print(raw_slice)  # This is a pyarrow.Table with the slice of data
    
    # Close the query after processing
    query.close()

