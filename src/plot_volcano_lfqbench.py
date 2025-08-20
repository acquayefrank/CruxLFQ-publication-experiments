import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# --- Data loading and merging (same as plot_fig2.py) ---
base = "../results/formatted_results/"
crux = base + "crux-lfq-mod-pep.txt_formatted"
flash = base + "FlashLFQ+mods+protein_id_modpep.txt_formatted"
ion = base + "ionquant_combined_modified_peptide.tsv_formatted"
maxq = base + "maxquant_peptides.txt_formatted"
sage = base + "sage_lfq.tsv_formatted"
proteomics = base + "proteomicslfq.mzTab_formatted"

df_crux = pd.read_csv(crux, sep="\t")
df_flash = pd.read_csv(flash, sep="\t")
df_sage = pd.read_csv(sage, sep="\t")
df_ion = pd.read_csv(ion, sep="\t")
df_max_maxquant = pd.read_csv(maxq, sep="\t")
df_proteomics = pd.read_csv(proteomics, sep="\t")

df_sage = df_sage.fillna(0)

df1 = pd.merge(df_crux, df_flash, left_on="id", right_on="id", how="inner", suffixes=("_crux", "_flash"))
df2 = pd.merge(df1, df_sage, left_on="id", right_on="id", how="inner", suffixes=("_df1", "_sage"))
df3 = pd.merge(df2, df_ion, left_on="id", right_on="id", how="inner", suffixes=("_df2", "_ion"))
df4 = pd.merge(df3, df_max_maxquant, left_on="id", right_on="id", how="inner", suffixes=("_df3", "_maxquant"))
df = pd.merge(df4, df_proteomics, left_on="id", right_on="id", how="inner", suffixes=("_df4", "_proteomics"))

# --- Column mappings (same as plot_fig2.py) ---
crux_cols_map = {
    "A1": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_06_161103_A1_HCD_OT_4ul.mzML",
    "A2": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_07_161103_A2_HCD_OT_4ul.mzML",
    "A3": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_16_161103_A3_HCD_OT_4ul.mzML",
    "A4": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_17_161103_A4_HCD_OT_4ul.mzML",
    "B1": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML",
    "B2": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_08_161103_B2_HCD_OT_4ul.mzML",
    "B3": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_15_161103_B3_HCD_OT_4ul.mzML",
    "B4": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_18_161103_B4_HCD_OT_4ul.mzML",
    "C1": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_24_161103_C1_HCD_OT_4ul.mzML",
    "C2": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_09_161103_C2_HCD_OT_4ul.mzML",
    "C3": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_14_161103_C3_HCD_OT_4ul.mzML",
    "C4": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_19_161103_C4_HCD_OT_4ul.mzML",
    "D1": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_22_161103_D1_HCD_OT_4ul.mzML",
    "D2": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_11_161103_D2_HCD_OT_4ul.mzML",
    "D3": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_12_161103_D3_HCD_OT_4ul.mzML",
    "D4": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_21_161103_D4_HCD_OT_4ul.mzML",
    "E1": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_23_161103_E1_HCD_OT_4ul.mzML",
    "E2": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_10_161103_E2_HCD_OT_4ul.mzML",
    "E3": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_13_161103_E3_HCD_OT_4ul.mzML",
    "E4": "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_20_161103_E4_HCD_OT_4ul.mzML",
}
flash_cols_map = {
    "A1": "Intensity_B02_06_161103_A1_HCD_OT_4ul",
    "A2": "Intensity_B02_07_161103_A2_HCD_OT_4ul",
    "A3": "Intensity_B02_16_161103_A3_HCD_OT_4ul",
    "A4": "Intensity_B02_17_161103_A4_HCD_OT_4ul",
    "B1": "Intensity_B02_001_161103_B1_HCD_OT_4ul",
    "B2": "Intensity_B02_08_161103_B2_HCD_OT_4ul",
    "B3": "Intensity_B02_15_161103_B3_HCD_OT_4ul",
    "B4": "Intensity_B02_18_161103_B4_HCD_OT_4ul",
    "C1": "Intensity_B02_24_161103_C1_HCD_OT_4ul",
    "C2": "Intensity_B02_09_161103_C2_HCD_OT_4ul",
    "C3": "Intensity_B02_14_161103_C3_HCD_OT_4ul",
    "C4": "Intensity_B02_19_161103_C4_HCD_OT_4ul",
    "D1": "Intensity_B02_22_161103_D1_HCD_OT_4ul",
    "D2": "Intensity_B02_11_161103_D2_HCD_OT_4ul",
    "D3": "Intensity_B02_12_161103_D3_HCD_OT_4ul",
    "D4": "Intensity_B02_21_161103_D4_HCD_OT_4ul",
    "E1": "Intensity_B02_23_161103_E1_HCD_OT_4ul",
    "E2": "Intensity_B02_10_161103_E2_HCD_OT_4ul",
    "E3": "Intensity_B02_13_161103_E3_HCD_OT_4ul",
    "E4": "Intensity_B02_20_161103_E4_HCD_OT_4ul",
}
sage_cols_map = {
    "A1": "B02_06_161103_A1_HCD_OT_4ul.mzML", "A2": "B02_07_161103_A2_HCD_OT_4ul.mzML", "A3": "B02_16_161103_A3_HCD_OT_4ul.mzML", "A4": "B02_17_161103_A4_HCD_OT_4ul.mzML",
    "B1": "B02_001_161103_B1_HCD_OT_4ul.mzML", "B2": "B02_08_161103_B2_HCD_OT_4ul.mzML", "B3": "B02_15_161103_B3_HCD_OT_4ul.mzML", "B4": "B02_18_161103_B4_HCD_OT_4ul.mzML",
    "C1": "B02_24_161103_C1_HCD_OT_4ul.mzML", "C2": "B02_09_161103_C2_HCD_OT_4ul.mzML", "C3": "B02_14_161103_C3_HCD_OT_4ul.mzML", "C4": "B02_19_161103_C4_HCD_OT_4ul.mzML",
    "D1": "B02_22_161103_D1_HCD_OT_4ul.mzML", "D2": "B02_11_161103_D2_HCD_OT_4ul.mzML", "D3": "B02_12_161103_D3_HCD_OT_4ul.mzML", "D4": "B02_21_161103_D4_HCD_OT_4ul.mzML",
    "E1": "B02_23_161103_E1_HCD_OT_4ul.mzML", "E2": "B02_10_161103_E2_HCD_OT_4ul.mzML", "E3": "B02_13_161103_E3_HCD_OT_4ul.mzML", "E4": "B02_20_161103_E4_HCD_OT_4ul.mzML",
}
max_cols_map = { # MaxLFQ
    "A1": "Exp1_1 MaxLFQ Intensity", "A2": "Exp1_2 MaxLFQ Intensity",
    "A3": "Exp1_3 MaxLFQ Intensity", "A4": "Exp1_4 MaxLFQ Intensity",
    "B1": "Exp2_1 MaxLFQ Intensity", "B2": "Exp2_2 MaxLFQ Intensity",
    "B3": "Exp2_3 MaxLFQ Intensity", "B4": "Exp2_4 MaxLFQ Intensity",
    "C1": "Exp3_1 MaxLFQ Intensity", "C2": "Exp3_2 MaxLFQ Intensity",
    "C3": "Exp3_3 MaxLFQ Intensity", "C4": "Exp3_4 MaxLFQ Intensity",
    "D1": "Exp4_1 MaxLFQ Intensity", "D2": "Exp4_2 MaxLFQ Intensity",
    "D3": "Exp4_3 MaxLFQ Intensity", "D4": "Exp4_4 MaxLFQ Intensity",
    "E1": "Exp5_1 MaxLFQ Intensity", "E2": "Exp5_2 MaxLFQ Intensity",
    "E3": "Exp5_3 MaxLFQ Intensity", "E4": "Exp5_4 MaxLFQ Intensity",
}
ion_cols_map = {
    "A1": "A_1 Intensity", "A2": "A_2 Intensity",
    "A3": "A_3 Intensity", "A4": "A_4 Intensity",
    "B1": "B_1 Intensity", "B2": "B_2 Intensity",
    "B3": "B_3 Intensity", "B4": "B_4 Intensity",
    "C1": "C_1 Intensity", "C2": "C_2 Intensity",
    "C3": "C_3 Intensity", "C4": "C_4 Intensity",
    "D1": "D_1 Intensity", "D2": "D_2 Intensity",
    "D3": "D_3 Intensity", "D4": "D_4 Intensity",
    "E1": "E_1 Intensity", "E2": "E_2 Intensity",
    "E3": "E_3 Intensity", "E4": "E_4 Intensity",
}
maxquant_cols_map = { # MaxQuant directly
    "A1": "Intensity A1", "A2": "Intensity A2",
    "A3": "Intensity A3", "A4": "Intensity A4",
    "B1": "Intensity B1", "B2": "Intensity B2",
    "B3": "Intensity B3", "B4": "Intensity B4",
    "C1": "Intensity C1", "C2": "Intensity C2",
    "C3": "Intensity C3", "C4": "Intensity C4",
    "D1": "Intensity D1", "D2": "Intensity D2",
    "D3": "Intensity D3", "D4": "Intensity D4",
    "E1": "Intensity E1", "E2": "Intensity E2",
    "E3": "Intensity E3", "E4": "Intensity E4",
}
proteomics_cols_map = {
    "A1": "A1", "A2": "A2", "A3": "A3", "A4": "A4",
    "B1": "B1", "B2": "B2", "B3": "B3", "B4": "B4",
    "C1": "C1", "C2": "C2", "C3": "C3", "C4": "C4",
    "D1": "D1", "D2": "D2", "D3": "D3", "D4": "D4",
    "E1": "E1", "E2": "E2", "E3": "E3", "E4": "E4",
}

# --- LFQbench configuration (matching R code) ---
LFQBENCH_CONFIG = {
    'PlotPointSize': 0.5,
    'LogRatioPlotRange': (-3, 3),
    'LogRatioValidityRangeSDFactor': 50,
    'SignificanceThreshold': 0.05,
    'SpeciesTags': {
        'HUMAN': '_HUMAN',
        'ECOLI': '_ECOLI'
    },
    'SampleComposition': {
        'species': ['HUMAN', 'ECOLI'],
        'A': [1, 1],
        'B': [1, 1.5], 
        'C': [1, 2],
        'D': [1, 2.5],
        'E': [1, 3]
    }
}

# # Filter for peptides quantified by all methods (no missing values)
# # Identify all intensity columns from each method
# all_intensity_cols = []
# all_intensity_cols += list(crux_cols_map.values())
# all_intensity_cols += list(flash_cols_map.values())
# all_intensity_cols += list(sage_cols_map.values())
# all_intensity_cols += list(max_cols_map.values())
# all_intensity_cols += list(ion_cols_map.values())
# all_intensity_cols += list(maxquant_cols_map.values())
# all_intensity_cols += list(proteomics_cols_map.values())
# # Only keep columns that exist in the merged DataFrame
# all_intensity_cols = [col for col in all_intensity_cols if col in df.columns]
# # Filter: keep only rows where all intensity columns are not null and >0
# df = df.dropna(subset=all_intensity_cols)
# df = df[(df[all_intensity_cols] > 0).all(axis=1)]
# ###############################################################

# --- LFQbench-style MA plot function (matching original plot axes) ---
def ma_plot_data_lfqbench(df, cols_map, groupA_keys, groupB_keys, group_name, min_valid_per_group=1):
    """
    LFQbench-style MA plot data generation matching R package approach
    X-axis: Log2(group B average intensity) - specific to each comparison
    Y-axis: Log2 fold change (A vs B)
    """
    groupA = [cols_map[k] for k in groupA_keys if cols_map[k] in df.columns]
    groupB = [cols_map[k] for k in groupB_keys if cols_map[k] in df.columns]
    
    print(f"GroupA columns found: {len(groupA)} out of {len(groupA_keys)}")
    print(f"GroupB columns found: {len(groupB)} out of {len(groupB_keys)}")
    
    if not groupA or not groupB:
        print("Warning: No matching columns found!")
        return np.array([]), np.array([]), 0, [], group_name
    
    # Replace 0 with NaN for proper handling of missing values
    dataA = df[groupA].replace(0, np.nan)
    dataB = df[groupB].replace(0, np.nan)
    
    # Create masks for different peptide types (simulating MBR vs direct identification)
    # Check how many samples have valid data for each peptide
    valid_countA = dataA.notna().sum(axis=1)
    valid_countB = dataB.notna().sum(axis=1)
    
    # Basic filtering: require at least 1 valid value per group
    maskA = valid_countA >= min_valid_per_group
    maskB = valid_countB >= min_valid_per_group
    mask = maskA & maskB
    
    print(f"Peptides passing filter: {mask.sum()} out of {len(df)}")
    
    if mask.sum() == 0:
        print("Warning: No peptides passed filtering!")
        return np.array([]), np.array([]), 0, [], group_name
    
    dataA_filtered = dataA[mask]
    dataB_filtered = dataB[mask]
    valid_countA_filtered = valid_countA[mask]
    valid_countB_filtered = valid_countB[mask]
    
    # Calculate means for each group (raw intensities)
    meanA = dataA_filtered.mean(axis=1, skipna=True)
    meanB = dataB_filtered.mean(axis=1, skipna=True)
    
    # Remove any remaining NaN or infinite values
    valid_mask = np.isfinite(meanA) & np.isfinite(meanB) & (meanA > 0) & (meanB > 0)
    meanA = meanA[valid_mask]
    meanB = meanB[valid_mask]
    valid_countA_final = valid_countA_filtered[valid_mask]
    valid_countB_final = valid_countB_filtered[valid_mask]
    
    if len(meanA) == 0:
        print("Warning: No valid intensities after filtering!")
        return np.array([]), np.array([]), 0, [], group_name
    
    # X-axis: Log2 of group B average intensity (specific to each comparison)
    log2_groupB_intensity = np.log2(meanB)
    
    # Y-axis: Log2 fold change (A vs B)
    log2_fold_change = np.log2(meanA) - np.log2(meanB)
    
    # Classify peptides by species based on protein IDs
    # Green: Human peptides, Orange: E.coli peptides
    df_filtered = df[mask]
    peptide_types = []
    for i in range(len(meanA)):
        idx = df_filtered.index[valid_mask][i]
        # Check if 'Protein ID' column exists, otherwise fall back to 'id'
        
        if 'Protein ID_df2' in df_filtered.columns:
            protein_id = str(df_filtered.loc[idx, 'Protein ID_df2']).upper()
        elif 'id' in df_filtered.columns:
            protein_id = str(df_filtered.loc[idx, 'id']).upper()
        else:
            protein_id = ""
        
        if '_HUMAN' in protein_id or 'HUMAN' in protein_id:
            peptide_types.append('human')
        elif '_ECOLI' in protein_id or 'ECOLI' in protein_id or 'COLI' in protein_id:
            peptide_types.append('ecoli')
        else:
            peptide_types.append('unknown')
    
    print(f"Final data points: {len(log2_groupB_intensity)}")
    print(f"Peptide species: {peptide_types.count('human')} human (green), {peptide_types.count('ecoli')} E.coli (orange), {peptide_types.count('unknown')} unknown")
    
    return log2_groupB_intensity, log2_fold_change, len(log2_groupB_intensity), peptide_types, group_name

# --- LFQbench-style MA plotting function ---
def create_lfqbench_ma_plots():
    """
    Create MA plots matching LFQbench R package style exactly:
    - X-axis: Log₁₀(average intensity)
    - Y-axis: Log₁₀ fold change
    - Methods as rows, Comparisons as columns
    - Green and orange colors matching original
    """
    
    # methods = [
    #     ('CruxLFQ', crux_cols_map),
    #     ('FlashLFQ', flash_cols_map),
    #     ('SageLFQ', sage_cols_map),
    #     ('IonQuant', ion_cols_map),
    #     ('MaxQuant', maxquant_cols_map),
    #     ('ProteomicsLFQ', proteomics_cols_map)
    # ]

    methods = [
        ('IonQuant', ion_cols_map),
        ('MaxQuant', maxquant_cols_map),
        ('FlashLFQ', flash_cols_map),
    ]
    
    group_pairs = [
        (['A1', 'A2', 'A3', 'A4'], ['B1', 'B2', 'B3', 'B4'], 'Human 1:1\nEcoli 1:1.5', 'B'),
        (['A1', 'A2', 'A3', 'A4'], ['C1', 'C2', 'C3', 'C4'], 'Human 1:1\nEcoli 1:2', 'C'),
        (['A1', 'A2', 'A3', 'A4'], ['D1', 'D2', 'D3', 'D4'], 'Human 1:1\nEcoli 1:2.5', 'D'),
        (['A1', 'A2', 'A3', 'A4'], ['E1', 'E2', 'E3', 'E4'], 'Human 1:1\nEcoli 1:3', 'E'),
    ]
    
    # Create figure with layout: 6 methods × 4 comparisons
    fig, axes = plt.subplots(len(methods), len(group_pairs), 
                            figsize=(20, 24), 
                            sharey=False, sharex=False)
    
    # Ensure axes is always 2D array
    if len(group_pairs) == 1:
        axes = axes[:, None]
    if len(methods) == 1:
        axes = axes[None, :]
    
    # Set style matching original
    plt.style.use('default')
    fig.patch.set_facecolor('white')
    
    for row, (method, cmap) in enumerate(methods):
        for col, (groupA, groupB, label, group_letter) in enumerate(group_pairs):
            ax = axes[row, col]
            
            # Get data using MA plot calculation
            print(f"\nProcessing {method} - {label}:")
            log2_groupB_intensity, log2_fold_change, n_peptides, peptide_types, group_name = ma_plot_data_lfqbench(df, cmap, groupA, groupB, group_letter)
            
            if len(log2_groupB_intensity) == 0:
                print(f"No data for {method} - {label}")
                ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes, 
                       ha='center', va='center', fontsize=12, color='red')
                # Still set up the axes even with no data
                ax.set_xlim(18, 30)
                ax.set_ylim(-3, 3)
                ax.grid(True, alpha=0.3, linewidth=0.5, color='lightgray')
                ax.set_facecolor('white')
                continue
            
            # Create scatter plot with species-based colors
            # Green for human peptides, Orange for E.coli peptides
            human_mask = np.array([t == 'human' for t in peptide_types])
            ecoli_mask = np.array([t == 'ecoli' for t in peptide_types])
            unknown_mask = np.array([t == 'unknown' for t in peptide_types])
            
            # Determine which dataset is larger and plot larger one first (bottom), smaller one on top
            human_count = human_mask.sum()
            ecoli_count = ecoli_mask.sum()
            unknown_count = unknown_mask.sum()
            
            # Create list of (count, mask, color, species_label) tuples and sort by count (descending)
            datasets = []
            if unknown_count > 0:
                datasets.append((unknown_count, unknown_mask, '#808080', 'Unknown'))
            if ecoli_count > 0:
                datasets.append((ecoli_count, ecoli_mask, '#FF8C00', 'E.coli'))
            if human_count > 0:
                datasets.append((human_count, human_mask, '#228B22', 'Human'))
            
            # Sort by count (largest first) so smaller datasets are plotted on top
            datasets.sort(key=lambda x: x[0], reverse=True)
            
            # Plot in order: largest to smallest (smaller datasets on top)
            for count, mask, color, species_label in datasets:
                ax.scatter(log2_groupB_intensity[mask], log2_fold_change[mask], 
                          alpha=0.7, s=2, c=color, edgecolors='none', label=f'{species_label} (n={count})')
            
            # Add horizontal reference line at y=0
            ax.axhline(y=0, color='gray', linestyle='-', lw=0.5, alpha=0.7)
            
            # Set axis limits matching the original image exactly
            ax.set_xlim(18, 30)  # X-axis range from original
            ax.set_ylim(-3, 3)   # Y-axis range from original
            
            # Grid and styling to match original exactly
            ax.grid(True, alpha=0.2, linewidth=0.5, color='lightgray')
            ax.set_facecolor('white')  # White background instead of light gray
            
            # Remove top and right spines to match original
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.8)
            ax.spines['bottom'].set_linewidth(0.8)
            ax.spines['left'].set_color('black')
            ax.spines['bottom'].set_color('black')
            
            # Labels and titles matching original exactly
            if row == 0:
                ax.set_title(label, fontsize=11, fontweight='bold', pad=15)
            
            # Method names on the left - moved further left to prevent overlap
            if col == 0:
                ax.text(-0.18, 0.5, method, 
                       transform=ax.transAxes, rotation=90, 
                       va='center', ha='center', fontsize=10, fontweight='bold')
            
            # Axis labels - specific to each comparison group
            if row == len(methods) - 1:
                ax.set_xlabel(f'Log₂({group_letter})', fontsize=10, fontweight='bold')
            else:
                ax.set_xticklabels([])
                
            if col == 0:
                ax.set_ylabel('Log₂ fold change', fontsize=9, fontweight='bold')
            else:
                ax.set_yticklabels([])
            
            # Tick parameters matching original
            ax.tick_params(axis='both', which='major', labelsize=9, 
                          length=4, width=0.8, color='black')
            
            # Add tick marks on all sides like original
            ax.tick_params(axis='x', which='major', top=True, bottom=True)
            ax.tick_params(axis='y', which='major', left=True, right=True)
            
            # Add number of points text in top-right corner
            if len(log2_groupB_intensity) > 0:
                ax.text(0.95, 0.95, f'n={n_peptides}', 
                       transform=ax.transAxes, fontsize=8, 
                       ha='right', va='top', 
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Main title matching original exactly
    fig.suptitle('Peptide-Level Quantification', fontsize=16, fontweight='bold', y=0.96)
    
    # Adjust layout to match original spacing
    plt.tight_layout(rect=[0.08, 0.03, 1, 0.93])  # Increased left margin for method labels
    plt.subplots_adjust(wspace=0.15, hspace=0.20)  # Increased hspace to prevent y-axis text overlap
    
    # Save plots
    plt.savefig('lfqbench_ma_plots_matching_original.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('lfqbench_ma_plots_matching_original.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print("LFQbench-style MA plots saved (matching original image exactly)")
    print(f"Layout: {len(methods)} methods × {len(group_pairs)} comparisons")
    print("X-axis: Log₂(B), Log₂(C), Log₂(D), Log₂(E) - specific to each comparison group")
    print("Y-axis: Log₂ fold change") 
    print("Colors: Green (human peptides), Orange (E.coli peptides)")
    print(f"Methods included: {', '.join([method[0] for method in methods])}")
    
    # Generate summary statistics
    print("\nSummary Statistics:")
    print("-" * 60)
    for method_name, cols_map in methods:
        print(f"\n{method_name}:")
        for groupA, groupB, label, group_letter in group_pairs:
            _, _, n_peptides, peptide_types, _ = ma_plot_data_lfqbench(df, cols_map, groupA, groupB, group_letter)
            if n_peptides > 0:
                human_count = sum(1 for t in peptide_types if t == 'human')
                ecoli_count = sum(1 for t in peptide_types if t == 'ecoli')
                unknown_count = sum(1 for t in peptide_types if t == 'unknown')
                print(f"  {group_letter}: {n_peptides} peptides ({human_count} human, {ecoli_count} E.coli, {unknown_count} unknown)")
            else:
                print(f"  {group_letter}: No data")

if __name__ == "__main__":
    print("Creating LFQbench-style MA plots matching original image...")
    print("X-axis: Log₁₀(average intensity)")
    print("Y-axis: Log₁₀ fold change")
    print("Methods: CruxLFQ, FlashLFQ, SageLFQ, IonQuant, MaxQuant, ProteomicsLFQ")
    print("=" * 80)
    
    # Print column availability for debugging
    print("\nDebugging column mappings:")
    print(f"DataFrame columns: {len(df.columns)} total")
    print(f"Sample CruxLFQ columns in df: {[col for col in crux_cols_map.values() if col in df.columns][:3]}...")
    print(f"Sample IonQuant columns in df: {[col for col in ion_cols_map.values() if col in df.columns][:3]}...")
    
    create_lfqbench_ma_plots()
    
    print("\nMA plot analysis complete!")
    print("Axes now match the original R LFQbench output")