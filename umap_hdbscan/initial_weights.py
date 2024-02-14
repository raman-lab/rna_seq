# Generates the reference dictionary which is used to caluclate fold enrichment of amino acids at each position of the variants after induction with the different ligands. DNA counts from each replkicate for each ligand was used. 
# 1. For each ligand, find files /FLC_{lignd}_rep{1-3}_barcode_all_counts.csv
# 2. Get "DNA_FLC_{ligand}_*} counts along with the corresponding barcode
# 3. Map each barcode back to a variant using /Volumes/sraman4/General/Nate/TtgR/data_processing/seq_data_full/FLC/sum1_withwt/FLC_Nal/funclib_all_barcodes_withwt.csv, summing up all the barcodes belonging to unique variants.
# 4. Map each variant back to amino acid sequence using /Volumes/Jackie/TtgR_RNASeq/umap/avg_foldenrichment_full.csv
# 5. Get the frequencies of each amino acid at each position using the weights determined by the DNA counts.

# Current working directory must be the folder containing the folders for each ligand.


import pandas as pd
import os, re
from collections import defaultdict




ligands = ligands = ['4Hy', 'End', 'Ndes', 'Tam', 'Phlo', 'Nar', 'Quin', 'Nal', 'EllA']

barcodemap = pd.read_csv("funclib_all_barcodes_withwt.csv")
bcmap_dict = {}
# Create a mask of non-NaN values in the DataFrame
non_nan_mask = barcodemap.notna()
# Iterate through the columns (excluding the first column)
for col_name in barcodemap.columns[1:]:
    # Create a mask for rows where the column is not NaN
    mask = non_nan_mask[col_name]
    
    # Add the non-NaN values as keys to the dictionary
    values = barcodemap.loc[mask, col_name].tolist()
    variants = barcodemap.loc[mask, "variants"].tolist()
    
    bcmap_dict.update(dict(zip(values, variants)))

# Variant to mutated positions mapping
aa_map = pd.read_csv("avg_foldenrichment_full.csv")
variant_aa_dict = dict(zip(aa_map['Variant_ID'], aa_map['Mut_Pos']))

filename_pattern = re.compile(r"FLC_(\w+)_rep[1-3]_barcode_all_counts.csv")
lig_ref_weights = {}
current_directory = os.getcwd() 

# Loop through ligands to find files containing DNA barcode counts
for ligand in ligands:
    dataframes = []
    columns = []
    variant_counts_ls = []
    folder_name = f"FLC_{ligand}"
    folder_path = os.path.join(current_directory, folder_name)
    
    #From thos files that contain the DNA barcode counts, extract the relevant column and add to a list
    for filename in os.listdir(folder_path):
        match = filename_pattern.match(filename)
        if match:
            file_ligand = match.group(1)
            if file_ligand == ligand:
                file_path = os.path.join(folder_path, filename)
                df = pd.read_csv(file_path, index_col=0)
                for col in df.columns:
                    if col.startswith(f'DNA_FLC_{ligand}'):
                        columns.append(df[col])
    
    # Loop through the list of relevant columns and get counts for the funclib variants and store in dictionary
    for series in columns:
        variant_counts = defaultdict(int)
        for barcode, count in series.items():
            try: # Skip barcodes which cannot be mapped
                variant_counts[variant_aa_dict[bcmap_dict[barcode]]] += count
            except:
                continue
        variant_counts_ls.append(variant_counts)
    ## Output a csv with the DNA counts for each of the variants
    avg_dna_count = defaultdict(list)
    for d in variant_counts_ls:
        for k, v in d.items():
            avg_dna_count[k].append(v)
    avg_dna_count_dict = {key: sum(v)/len(v) for key, v in avg_dna_count.items()}
    avg_dna_count_df = pd.DataFrame(avg_dna_count_dict.items(), columns=['Mut_Pos', 'DNA_Counts'])
    avg_dna_count_df.to_csv('dna_counts_full.csv')
        
    # Gather positional counts in a new dictionary for each replicate
    mut_count_ls = []
    for variant_counts in variant_counts_ls:
        mut_count = {}
        for mutpos, count in variant_counts.items():
        # Iterate through each position in the string
            for position, aa in enumerate(mutpos):
            # Create an inner dictionary for the position if it doesn't exist
                if position+1 not in mut_count:
                    mut_count[position+1] = {}
            # Update the inner dictionary with the character and its summed value
                if aa in mut_count[position+1]:
                    mut_count[position+1][aa] += count
                else:
                    mut_count[position+1][aa] = count
        mut_count_ls.append(mut_count)
    
    # Ger the frequencies by replicate 
    for mut_count in mut_count_ls:
        for pos, inner_dict in mut_count.items():
            total = sum(inner_dict.values())
            for aa, count in inner_dict.items():
                inner_dict[aa] = count/total 
                
    # Create a new dicitoanry using the averages from the three replicates
    averages_dict = {}
    # Iterate through each dictionary in the list
    for d in mut_count_ls:
        # Iterate through each key in the nested dictionary
        for key, inner_dict in d.items():
            if key not in averages_dict:
                averages_dict[key] = {}  # Create a new inner dictionary
            
            # Iterate through each inner key and value
            for inner_key, value in inner_dict.items():
                if inner_key not in averages_dict[key]:
                    averages_dict[key][inner_key] = []  # Create a list to store values
                
                averages_dict[key][inner_key].append(value)  # Store the value in the list
    
    # Calculate the average for each inner key in each key
    for key, inner_dict in averages_dict.items():
        for inner_key, values in inner_dict.items():
            averages_dict[key][inner_key] = sum(values) / len(values)   
    
    # Sort the inner dict by alphabetical order
    for key, inner_dict in averages_dict.items():
        averages_dict[key] = dict(sorted(inner_dict.items()))
    lig_ref_weights[ligand] = averages_dict
