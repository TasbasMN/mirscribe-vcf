import os
import pandas as pd
from tqdm import tqdm

df = pd.read_csv("data/560.txt", sep="\t")
cols_to_drop = ["Project", "ID", "Genome", "mut_type", "Type"]
df = df.drop(columns=cols_to_drop)

# setting vcfs to exclude
target_folder = "data/sample_vcfs"
vcfs_to_exclude = {os.path.basename(x).split(".")[0] for x in os.listdir(target_folder) if x.startswith("PD")}

output_dir = "data/560"
os.makedirs(output_dir, exist_ok=True)


for name, group in tqdm(df.groupby("Sample"), total=len(df["Sample"].unique())):
    if name in vcfs_to_exclude:
        continue
    
    # Create the appropriate directory
    if name.startswith(tuple(f"PD{i}" for i in range(1, 10))):
        sample_dir = os.path.join(output_dir, name[:3])
        os.makedirs(sample_dir, exist_ok=True)
    else:
        sample_dir = output_dir
    
    # Save the VCF file
    file_path = os.path.join(sample_dir, f"{name}.vcf")
    group.to_csv(file_path, sep="\t", index=False, header=False)









