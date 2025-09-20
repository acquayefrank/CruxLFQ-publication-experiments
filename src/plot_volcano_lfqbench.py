import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# --- Data loading and merging (same as plot_fig2.py) ---
base = "../results/formatted_results/"
crux = base + "crux-lfq-mod-pep-in.txt_formatted"
flash = base + "FlashLFQ+mods+protein_id_modpep-in.txt_formatted"
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

# --- Expected fold changes for each comparison ---
def get_expected_fold_change(group_letter):
    """
    Get the expected log2 fold change for E.coli peptides based on the experimental design
    Human: always 1:1 (no change) = 0 on log2 scale
    E.coli: varies by group
    """
    expected_ratios = {
        'B': 1.5,   # E.coli 1.5:1 = log2(1/1.5) = -0.585
        'C': 2.0,   # E.coli 2:1 = log2(1/2) = -1.0
        'D': 2.5,   # E.coli 2.5:1 = log2(1/2.5) = -1.322
        'E': 3.0    # E.coli 3:1 = log2(1/3) = -1.585
    }
    return np.log2(1 / expected_ratios.get(group_letter, 1.0))

# --- Enhanced LFQbench-style MA plotting function with box plots ---
def create_lfqbench_ma_plots_with_boxplots():
    """
    Create MA plots with box plots on the right side
    Left: MA plot (X: Log₂ intensity, Y: Log₂ fold change)
    Right: Box plot showing fold change distributions by species
    """
    
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
    
    # Create figure with subplots: each "cell" has 2 subplots (MA plot + box plot)
    fig = plt.figure(figsize=(28, 18))
    
    # Create a grid: 3 rows (methods) x 8 columns (4 comparisons × 2 subplots each)
    n_methods = len(methods)
    n_comparisons = len(group_pairs)
    
    # Use GridSpec for more control over subplot spacing
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(n_methods, n_comparisons * 2, 
                  figure=fig,
                  width_ratios=[3, 1] * n_comparisons,  # MA plot wider than box plot
                  hspace=0.25, wspace=0.15)
    
    fig.patch.set_facecolor('white')
    
    for row, (method, cmap) in enumerate(methods):
        for col, (groupA, groupB, label, group_letter) in enumerate(group_pairs):
            
            # Create MA plot (left subplot)
            ax_ma = fig.add_subplot(gs[row, col*2])
            # Create box plot (right subplot)  
            ax_box = fig.add_subplot(gs[row, col*2 + 1])
            
            # Get data using MA plot calculation
            print(f"\nProcessing {method} - {label}:")
            log2_groupB_intensity, log2_fold_change, n_peptides, peptide_types, group_name = ma_plot_data_lfqbench(df, cmap, groupA, groupB, group_letter)
            
            if len(log2_groupB_intensity) == 0:
                print(f"No data for {method} - {label}")
                # Handle no data case for both plots
                for ax in [ax_ma, ax_box]:
                    ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes, 
                           ha='center', va='center', fontsize=12, color='red')
                
                # Set up MA plot axes
                ax_ma.set_xlim(22, 36)
                ax_ma.set_ylim(-3, 3)
                ax_ma.grid(True, alpha=0.3, linewidth=0.5, color='lightgray')
                ax_ma.set_facecolor('white')
                
                # Set up box plot axes
                ax_box.set_xlim(-0.5, 1.5)
                ax_box.set_ylim(-3, 3)
                ax_box.grid(True, alpha=0.3, linewidth=0.5, color='lightgray')
                ax_box.set_facecolor('white')
                continue
            
            # ===== MA PLOT (LEFT) =====
            # Create scatter plot with species-based colors
            human_mask = np.array([t == 'human' for t in peptide_types])
            ecoli_mask = np.array([t == 'ecoli' for t in peptide_types])
            unknown_mask = np.array([t == 'unknown' for t in peptide_types])
            
            # Count peptides by species
            human_count = human_mask.sum()
            ecoli_count = ecoli_mask.sum()
            unknown_count = unknown_mask.sum()
            
            # Create list of (count, mask, color, species_label) and plot largest first
            datasets = []
            if unknown_count > 0:
                datasets.append((unknown_count, unknown_mask, '#808080', 'Unknown'))
            if ecoli_count > 0:
                datasets.append((ecoli_count, ecoli_mask, '#FF8C00', 'E.coli'))
            if human_count > 0:
                datasets.append((human_count, human_mask, '#228B22', 'Human'))
            
            # Sort by count (largest first) so smaller datasets are plotted on top
            datasets.sort(key=lambda x: x[0], reverse=True)
            
            # Plot scatter points
            for count, mask, color, species_label in datasets:
                ax_ma.scatter(log2_groupB_intensity[mask], log2_fold_change[mask], 
                            alpha=0.7, s=2, c=color, edgecolors='none', 
                            label=f'{species_label} (n={count})')
            
            # Add horizontal reference line at y=0
            ax_ma.axhline(y=0, color='gray', linestyle='-', lw=0.5, alpha=0.7)
            
            # Set MA plot limits and styling
            ax_ma.set_xlim(22, 36)
            ax_ma.set_xticks([22, 24, 26, 28, 30, 32, 34, 36])
            ax_ma.set_ylim(-3, 3)
            
            # Grid and styling
            ax_ma.grid(True, alpha=0.2, linewidth=0.5, color='lightgray')
            ax_ma.set_facecolor('white')
            ax_ma.spines['top'].set_visible(False)
            ax_ma.spines['right'].set_visible(False)
            ax_ma.spines['left'].set_linewidth(0.8)
            ax_ma.spines['bottom'].set_linewidth(0.8)
            
            # ===== BOX PLOT (RIGHT) =====
            # Prepare data for box plot
            box_data = []
            box_labels = []
            box_colors = []
            
            # Get expected fold change for this comparison
            expected_ecoli_fc = get_expected_fold_change(group_letter)
            
            if human_count > 0:
                box_data.append(log2_fold_change[human_mask])
                box_labels.append(f'Human\n(n={human_count})')
                box_colors.append('#228B22')
            
            if ecoli_count > 0:
                box_data.append(log2_fold_change[ecoli_mask])
                box_labels.append(f'E.coli\n(n={ecoli_count})')
                box_colors.append('#FF8C00')
            
            if len(box_data) > 0:
                # Create box plot
                bp = ax_box.boxplot(box_data, labels=box_labels, patch_artist=True,
                                    widths=0.7,  # Increase this value for wider boxes
                                   boxprops=dict(linewidth=1.5),
                                   medianprops=dict(linewidth=2, color='black'),
                                   whiskerprops=dict(linewidth=1.5),
                                   capprops=dict(linewidth=1.5),
                                   flierprops=dict(marker='o', markersize=3, alpha=0.6))
                
                # Color the boxes
                for patch, color in zip(bp['boxes'], box_colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
            
            # Add reference lines on box plot
            ax_box.axhline(y=0, color='gray', linestyle='-', lw=1, alpha=0.7, label='No change')
            if group_letter in ['B', 'C', 'D', 'E']:
                ax_box.axhline(y=expected_ecoli_fc, color='red', linestyle='--', lw=1, alpha=0.8,
                              label=f'Expected E.coli: {expected_ecoli_fc:.2f}')
            
            # Box plot styling
            ax_box.set_ylim(-3, 3)
            ax_box.grid(True, alpha=0.2, linewidth=0.5, color='lightgray')
            ax_box.set_facecolor('white')
            ax_box.spines['top'].set_visible(False)
            ax_box.spines['right'].set_visible(False)
            ax_box.spines['left'].set_linewidth(0.8)
            ax_box.spines['bottom'].set_linewidth(0.8)
            
            # Rotate box plot labels
            ax_box.tick_params(axis='x', rotation=45, labelsize=8)
            ax_box.tick_params(axis='y', labelsize=9)
            
            # ===== LABELS AND TITLES =====
            # Titles (only on top row)
            if row == 0:
                ax_ma.set_title(label, fontsize=11, fontweight='bold', pad=15)
            
            # Method names on the left
            if col == 0:
                ax_ma.text(-0.18, 0.5, method, 
                          transform=ax_ma.transAxes, rotation=90, 
                          va='center', ha='center', fontsize=10, fontweight='bold')
            
            # X-axis labels (only on bottom row)
            if row == len(methods) - 1:
                ax_ma.set_xlabel(f'Log₂({group_letter})', fontsize=10, fontweight='bold')
                ax_box.set_xlabel('Species', fontsize=9, fontweight='bold')
            else:
                ax_ma.set_xticklabels([])
                ax_box.set_xticklabels([])
            
            # Y-axis labels (only on leftmost column)
            if col == 0:
                ax_ma.set_ylabel('Log₂ fold change', fontsize=9, fontweight='bold')
            else:
                ax_ma.set_yticklabels([])
                ax_box.set_yticklabels([])
            
            # Tick parameters
            ax_ma.tick_params(axis='both', which='major', labelsize=9, 
                             length=4, width=0.8, color='black')
            ax_box.tick_params(axis='both', which='major', labelsize=8,
                              length=4, width=0.8, color='black')
            
            # Add tick marks on all sides
            ax_ma.tick_params(axis='x', which='major', top=True, bottom=True)
            ax_ma.tick_params(axis='y', which='major', left=True, right=True)
            
            # Add number of points text on MA plot
            if len(log2_groupB_intensity) > 0:
                ax_ma.text(0.95, 0.95, f'n={n_peptides}', 
                          transform=ax_ma.transAxes, fontsize=8, 
                          ha='right', va='top', 
                          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                   alpha=0.8, edgecolor='gray'))
    
    # Main title
    fig.suptitle('Peptide-Level Quantification with Distribution Analysis', 
                 fontsize=16, fontweight='bold', y=0.96)
    
    # Adjust layout
    plt.tight_layout(rect=[0.06, 0.03, 1, 0.93])
    
    # Save plots
    plt.savefig('lfqbench_ma_plots_with_boxplots.pdf', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('lfqbench_ma_plots_with_boxplots.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print("LFQbench-style MA plots with box plots saved")
    print(f"Layout: {len(methods)} methods × {len(group_pairs)} comparisons (MA + Box plots)")
    print("Left panels: MA plots (X: Log₂ intensity, Y: Log₂ fold change)")
    print("Right panels: Box plots showing fold change distributions by species")
    print("Colors: Green (human peptides), Orange (E.coli peptides)")
    
    # Generate detailed summary statistics
    print("\nDetailed Summary Statistics:")
    print("=" * 80)
    for method_name, cols_map in methods:
        print(f"\n{method_name}:")
        print("-" * 40)
        for groupA, groupB, label, group_letter in group_pairs:
            log2_groupB_intensity, log2_fold_change, n_peptides, peptide_types, _ = ma_plot_data_lfqbench(df, cols_map, groupA, groupB, group_letter)
            if n_peptides > 0:
                human_mask = np.array([t == 'human' for t in peptide_types])
                ecoli_mask = np.array([t == 'ecoli' for t in peptide_types])
                unknown_mask = np.array([t == 'unknown' for t in peptide_types])
                
                human_count = human_mask.sum()
                ecoli_count = ecoli_mask.sum()
                unknown_count = unknown_mask.sum()
                
                expected_ecoli_fc = get_expected_fold_change(group_letter)
                
                print(f"  {group_letter} vs A: {n_peptides} total peptides")
                print(f"    Human: {human_count} peptides")
                if human_count > 0:
                    human_fc = log2_fold_change[human_mask]
                    print(f"      Mean fold change: {human_fc.mean():.3f} ± {human_fc.std():.3f}")
                    print(f"      Median fold change: {np.median(human_fc):.3f}")
                
                print(f"    E.coli: {ecoli_count} peptides")
                if ecoli_count > 0:
                    ecoli_fc = log2_fold_change[ecoli_mask]
                    print(f"      Mean fold change: {ecoli_fc.mean():.3f} ± {ecoli_fc.std():.3f}")
                    print(f"      Median fold change: {np.median(ecoli_fc):.3f}")
                    print(f"      Expected fold change: {expected_ecoli_fc:.3f}")
                    print(f"      Difference from expected: {ecoli_fc.mean() - expected_ecoli_fc:.3f}")
                
                if unknown_count > 0:
                    print(f"    Unknown: {unknown_count} peptides")
            else:
                print(f"  {group_letter}: No data")

if __name__ == "__main__":
    print("Creating LFQbench-style MA plots with box plots...")
    print("Left panels: MA plots (X: Log₂ intensity, Y: Log₂ fold change)")
    print("Right panels: Box plots showing fold change distributions")
    print("Methods: IonQuant, MaxQuant, FlashLFQ")
    print("=" * 80)
    
    # Print column availability for debugging
    print("\nDebugging column mappings:")
    print(f"DataFrame columns: {len(df.columns)} total")
    print(f"Sample IonQuant columns in df: {[col for col in ion_cols_map.values() if col in df.columns][:3]}...")
    print(f"Sample MaxQuant columns in df: {[col for col in maxquant_cols_map.values() if col in df.columns][:3]}...")
    
    create_lfqbench_ma_plots_with_boxplots()
    
    print("\nMA plot with box plot analysis complete!")
    print("Box plots show the distribution of log2 fold changes for each species")
    print("Red dashed lines indicate expected E.coli fold changes based on experimental design")