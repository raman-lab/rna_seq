import matplotlib.pyplot as plt
import pandas as pd
import statistics as stats
import statsmodels.api as sm
import math

# Code used to identigy DNA-bound variants from sequencing

bcmap = pd.read_csv("funclib_all_barcodes_withwt.csv")
bcmap_dict = {}
# Create a mask of non-NaN values in the DataFrame
non_nan_mask = bcmap.notna()
# Iterate through the columns (excluding the first column)
for col_name in bcmap.columns[1:]:
    # Create a mask for rows where the column is not NaN
    mask = non_nan_mask[col_name]
    
    # Add the non-NaN values as keys to the dictionary
    values = bcmap.loc[mask, col_name].tolist()
    variants = bcmap.loc[mask, "variants"].tolist()
    
    bcmap_dict.update(dict(zip(values, variants)))

dmso1 = pd.read_csv("FLC_DMSO_rep1_barcode_all_counts.csv")
dmso2 = pd.read_csv("FLC_DMSO_rep2_barcode_all_counts.csv")
dmso3 = pd.read_csv("FLC_DMSO_rep3_barcode_all_counts.csv")

etoh1 = pd.read_csv("FLC_EtOH_rep1_barcode_all_counts.csv")
etoh2 = pd.read_csv("FLC_EtOH_rep2_barcode_all_counts.csv")
etoh3 = pd.read_csv("FLC_EtOH_rep3_barcode_all_counts.csv")

dmso_dict = {}
water_dict = {}
etoh_dict = {}

for df in [dmso1, dmso2, dmso3]:
    df = df.dropna(subset=[df.columns[3], df.columns[5]]) # Remove rows which have no DNA Counts
    df = df.fillna(0) # Fill in NaNs in the RNA columns with 0 
    df['Variant'] = df['0'].map(bcmap_dict) # Map barcodes to variants
    df = df.dropna(subset=['Variant']) # Remove rows that are NaN
    unique_variants = df['Variant'].unique() # Get unique variant names
    for variant in unique_variants: # Find unique variant names
        subset = df.loc[df['Variant'] == variant]
        dmso_fi = subset[df.columns[4]].sum()/subset[df.columns[5]].sum() # For each variant, get ratio of RNA sum and DNA sum for DMSO
        water_fi = subset[df.columns[2]].sum()/subset[df.columns[3]].sum() # For each variant, get ratio of RNA sum and DNA sum for Water
        if variant in dmso_dict:    # Add a key for each variant whose values are lists containing fold induction from each replicate
            dmso_dict[variant].append(dmso_fi)
        else:
            dmso_dict[variant] = [dmso_fi]
        if variant in water_dict:
            water_dict[variant].append(water_fi)
        else:
            water_dict[variant] = [water_fi]

for df in [etoh1, etoh2, etoh3]:
    df = df.dropna(subset=[df.columns[3], df.columns[5]]) # Remove rows which have no DNA Counts
    df = df.fillna(0) # Fill in NaNs in the RNA columns with 0 
    df['Variant'] = df['0'].map(bcmap_dict) # Map barcodes to variants
    df = df.dropna(subset=['Variant']) # Remove rows that are NaN
    unique_variants = df['Variant'].unique() # Get unique variant names
    for variant in unique_variants: # Find unique variant names
        subset = df.loc[df['Variant'] == variant]
        etoh_fi = subset[df.columns[4]].sum()/subset[df.columns[5]].sum() # For each variant, get ratio of RNA sum and DNA sum for EtOH
        if variant in etoh_dict:    # Add a key for each variant whose values are lists containing fold induction from each replicate
            etoh_dict[variant].append(etoh_fi)
        else:
            etoh_dict[variant] = [etoh_fi]

dmso_ci_dict = {key: [stats.mean(values)] + list(sm.stats.DescrStatsW(values).tconfint_mean(alpha=0.05)) for key, values in dmso_dict.items()}
water_ci_dict = {key: [stats.mean(values)] + list(sm.stats.DescrStatsW(values).tconfint_mean(alpha=0.05)) for key, values in water_dict.items()}
etoh_ci_dict = {key: [stats.mean(values)] + list(sm.stats.DescrStatsW(values).tconfint_mean(alpha=0.05)) for key, values in etoh_dict.items()}

fig, axes = plt.subplots(3, sharex=True)
solvents = ["DMSO", "Water", "EtOH"]
for i, ci_dict in enumerate([dmso_ci_dict, water_ci_dict, etoh_ci_dict]):
    wt_lower = ci_dict['WT'][1]
    wt_upper = ci_dict['WT'][2]
    count=0
    for key, values in ci_dict.items():
        if len(values) == 3:
        # Check if the second value in the list is greater than the third "WT" value
            if not math.isnan(values[1]) and values[1] > wt_upper:
                count += 1
            
    means = [values[0] for values in ci_dict.values()]
    axes[i].hist(means, bins=100, color='gray')
    axes[i].set_yscale('log')
    axes[i].axvline(wt_lower, color='red', linestyle='-')
    axes[i].axvline(wt_upper, color='red', linestyle='-')
    axes[i].set_title(solvents[i])
    
    axes[i].text(0.49, 0.85, f'Lower WT Bound: {round(wt_lower, 3)}', fontsize=12, transform=axes[i].transAxes)
    axes[i].text(0.49, 0.69, f'Upper WT Bound: {round(wt_upper, 3)}', fontsize=12, transform=axes[i].transAxes)
    axes[i].text(0.49, 0.53, f'Weak Repressor Count: {count}/{len(ci_dict)}', fontsize=12, transform=axes[i].transAxes)

plt.ylabel("Count")
plt.xlabel("RNA/DNA")
plt.tight_layout()
plt.savefig('fold_induction_histogram_95ci.pdf')
plt.show()
    
