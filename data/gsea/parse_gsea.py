import json
import pandas as pd

# Load the JSON data
with open("data/gsea/h.all.v2023.2.Hs.json") as file:
    gsea_data = json.load(file)

# Create a dictionary to store the data
data_dict = {}

# Iterate over each hallmark property
for hallmark, hallmark_data in gsea_data.items():
    # Extract the gene symbols for the current hallmark
    gene_symbols = hallmark_data['geneSymbols']
    
    # Add the gene symbols as a column in the data dictionary
    data_dict[hallmark] = pd.Series(True, index=gene_symbols)

# Create the DataFrame from the data dictionary
df = pd.DataFrame(data_dict)

# Fill missing values with False
df.fillna(False, inplace=True)

# Reset the index to make gene IDs as rows
df.reset_index(inplace=True)
df.rename(columns={'index': 'Gene ID'}, inplace=True)

df.to_csv('data/gsea/gsea.csv', index=False)