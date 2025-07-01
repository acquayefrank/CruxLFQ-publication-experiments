import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths (adjust as needed)
base = "../results/formatted_results/"
crux = base + "crux-lfq-mod-pep.txt_formatted"
flash = base + "FlashLFQ+mods+protein_id_modpep.txt_formatted"
ion = base + "ionquant_combined_peptide.tsv_formatted"
maxq = base + "maxquant_peptides.txt_formatted"
sage = base + "sage_lfq.tsv_formatted"
proteomics = base + "proteomicslfq.mzTab_formatted"

df_crux = pd.read_csv(crux, sep="\t")
df_flash = pd.read_csv(flash, sep="\t")
df_sage = pd.read_csv(sage, sep="\t", low_memory=False)
df_ion = pd.read_csv(ion, sep="\t")
df_max_maxquant = pd.read_csv(maxq, sep="\t")
df_proteomics = pd.read_csv(proteomics, sep="\t")

# Replace NaN values with 0 in df_sage
df_sage = df_sage.fillna(0)

df1 = pd.merge(
    df_crux,
    df_flash,
    left_on="id",
    right_on="id",
    how="inner",
    suffixes=("_crux", "_flash"),
)

df2 = pd.merge(
    df1,
    df_sage,
    left_on="id",
    right_on="id",
    how="inner",
    suffixes=("_df1", "_sage"),
)

df3 = pd.merge(
    df2,
    df_ion,
    left_on="id",
    right_on="id",
    how="inner",
    suffixes=("_df2", "_ion"),
)

df4 = pd.merge(
    df3,
    df_max_maxquant,
    left_on="id",
    right_on="id",
    how="inner",
    suffixes=("_df3", "_maxquant"),
)

df = pd.merge(
    df4,
    df_proteomics,
    left_on="id",
    right_on="id",
    how="inner",
    suffixes=("_df4", "_proteomics"),
)

# --- Category initialization and processing functions (unchanged) ---
def initialize_categories():
    categories = {}
    categories_list = [
        "A1", "A2", "A3", "A4",
        "B1", "B2", "B3", "B4",
        "C1", "C2", "C3", "C4",
        "D1", "D2", "D3", "D4",
        "E1", "E2", "E3", "E4",
    ]
    for category in categories_list:
        categories[category] = []
    return categories

def safe_divide(arr1, arr2):
    safe_array = [(ar1, ar2) for ar1, ar2 in zip(arr1, arr2) if ar1 != 0 and ar2 != 0]
    arr1_safe = [ar[0] for ar in safe_array]
    arr2_safe = [ar[1] for ar in safe_array]
    result = np.divide(arr1_safe, arr2_safe)
    return result

def calculate_data(a, b, c, d, e):
    return (
        safe_divide(np.array(a), np.array(b)),
        safe_divide(np.array(a), np.array(c)),
        safe_divide(np.array(a), np.array(d)),
        safe_divide(np.array(a), np.array(e)),
    )

def generate_categories(categories):
    a, b, c, d, e = [], [], [], [], []
    for key in categories:
        if "A" in key:
            a.extend(categories[key])
        elif "B" in key:
            b.extend(categories[key])
        elif "C" in key:
            c.extend(categories[key])
        elif "D" in key:
            d.extend(categories[key])
        elif "E" in key:
            e.extend(categories[key])
    return calculate_data(a, b, c, d, e)

# --- Optimized Data Splitting and Processing ---


# Pre-define column mappings for efficiency
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
    "A1": "Exp1_1 Intensity", "A2": "Exp1_2 Intensity",
    "A3": "Exp1_3 Intensity", "A4": "Exp1_4 Intensity",
    "B1": "Exp2_1 Intensity", "B2": "Exp2_2 Intensity",
    "B3": "Exp2_3 Intensity", "B4": "Exp2_4 Intensity",
    "C1": "Exp3_1 Intensity", "C2": "Exp3_2 Intensity",
    "C3": "Exp3_3 Intensity", "C4": "Exp3_4 Intensity",
    "D1": "Exp4_1 Intensity", "D2": "Exp4_2 Intensity",
    "D3": "Exp4_3 Intensity", "D4": "Exp4_4 Intensity",
    "E1": "Exp5_1 Intensity", "E2": "Exp5_2 Intensity",
    "E3": "Exp5_3 Intensity", "E4": "Exp5_4 Intensity",
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

# Create lists of column names for easier slicing
crux_cols = list(crux_cols_map.values())
flash_cols = list(flash_cols_map.values())
sage_cols = list(sage_cols_map.values())
max_cols = list(max_cols_map.values()) # MaxLFQ
ion_cols = list(ion_cols_map.values())
maxquant_cols = list(maxquant_cols_map.values()) # MaxQuant directly
proteomics_cols = list(proteomics_cols_map.values())

# Split the DataFrame once based on 'Protein ID_df2'
df_human = df[df["Protein ID_df2"].str.contains("HUMAN")].copy()
df_ecoli = df[df["Protein ID_df2"].str.contains("ECOLI")].copy()

def populate_categories_optimized(df_species, cols_map):
    categories = initialize_categories()
    for category_key, col_name in cols_map.items():
        if col_name in df_species.columns: # Check if column exists after merges
            categories[category_key].extend(df_species[col_name].tolist())
    return categories

human_crux_categories = populate_categories_optimized(df_human, crux_cols_map)
human_flash_categories = populate_categories_optimized(df_human, flash_cols_map)
human_sage_categories = populate_categories_optimized(df_human, sage_cols_map)
human_max_categories = populate_categories_optimized(df_human, max_cols_map)
human_ion_categories = populate_categories_optimized(df_human, ion_cols_map)
human_max_maxquant = populate_categories_optimized(df_human, maxquant_cols_map)
human_proteomics_categories = populate_categories_optimized(df_human, proteomics_cols_map)

ecoli_crux_categories = populate_categories_optimized(df_ecoli, crux_cols_map)
ecoli_flash_categories = populate_categories_optimized(df_ecoli, flash_cols_map)
ecoli_sage_categories = populate_categories_optimized(df_ecoli, sage_cols_map)
ecoli_max_categories = populate_categories_optimized(df_ecoli, max_cols_map)
ecoli_ion_categories = populate_categories_optimized(df_ecoli, ion_cols_map)
ecoli_max_maxquant = populate_categories_optimized(df_ecoli, maxquant_cols_map)
ecoli_proteomics_categories = populate_categories_optimized(df_ecoli, proteomics_cols_map)


# --- Generate data ---
HUMAN_CRUX = generate_categories(human_crux_categories)
HUMAN_FLASH = generate_categories(human_flash_categories)
HUMAN_SAGE = generate_categories(human_sage_categories)
HUMAN_MAX = generate_categories(human_max_categories)
HUMAN_ION = generate_categories(human_ion_categories)
HUMAN_MAX_MAXQUANT = generate_categories(human_max_maxquant)
HUMAN_PROTEOMICS = generate_categories(human_proteomics_categories)

ECOLI_CRUX = generate_categories(ecoli_crux_categories)
ECOLI_FLASH = generate_categories(ecoli_flash_categories)
ECOLI_SAGE = generate_categories(ecoli_sage_categories)
ECOLI_MAX = generate_categories(ecoli_max_categories)
ECOLI_ION = generate_categories(ecoli_ion_categories)
ECOLI_MAX_MAXQUANT = generate_categories(ecoli_max_maxquant)
ECOLI_PROTEOMICS = generate_categories(ecoli_proteomics_categories)

names = ["1.5-fold (B)", "2.0-fold (C)", "2.5-fold (D)", "3.0-fold (E)"]

# Helper to filter out NaNs from each array in a list
def filter_nan(arrays):
    return [np.array([x for x in arr if not np.isnan(x)]) for arr in arrays]

# --- Plotting functions (unchanged, as they are not the bottleneck) ---

# Filter out NaNs and log2 transform
def preprocess(data_tuple):
    return [np.array([np.log2(x) for x in l if x > 0 and not np.isnan(x)]) for l in data_tuple]

def plot_combined_boxplot(
    method_data_by_group, title, filename, y_label="log2 fold change"
):
    fig, ax = plt.subplots(figsize=(25, 17))

    num_methods = len(method_data_by_group)
    num_groups = len(names)
    box_width = 0.7
    spacing = 1.5  # space between boxes in a group
    group_spacing = spacing * num_methods + 1

    positions = []
    group_centers = []

    for group_idx in range(num_groups):
        start = group_idx * group_spacing
        pos = [start + i * spacing for i in range(num_methods)]
        positions.append(pos)
        group_centers.append(np.mean(pos))

    # Background shading
    custom_colors = ["#FFFFFF", "#D3D3D3", "#FFFFFF", "#D3D3D3"]
    for i, color in enumerate(custom_colors):
        ax.axvspan(
            group_centers[i] - spacing * 3,
            group_centers[i] + spacing * 3,
            facecolor=color,
            alpha=0.3,
        )

    # Plot each method
    colors = [
        "lightblue", "lightgreen", "lightyellow", "pink", "lightcoral", "lightgray"
    ]
    for method_idx in range(num_methods):
        data = [method_data_by_group[method_idx][i] for i in range(num_groups)]
        method_positions = [positions[i][method_idx] for i in range(num_groups)]
        ax.boxplot(
            data,
            positions=method_positions,
            widths=box_width,
            patch_artist=True,
            boxprops=dict(facecolor=colors[method_idx]),
            showfliers=False,
            showmeans=True,
        )

    # Fold-change reference line
    ax.axhline(y=0, color="b", linestyle="--", linewidth=0.7)

    # Add FC horizontal lines for ECOLI only
    if "ecoli" in filename.lower():
        x_min, _ = ax.get_xlim()
        # label_x = group_centers[0] - spacing * 5  # inside the first group
        # Dynamically adjust label_x to be slightly to the left of the first group's center
        label_x = min(group_centers) - (group_centers[1] - group_centers[0]) / 2 - 1.5


        fold_changes = [1.5, 2.0, 2.5, 3.0]
        for fc in fold_changes:
            log2_fc = np.log2(1 / fc)
            ax.axhline(y=log2_fc, color="r", linestyle="--", linewidth=0.7)
            ax.text(
                x=max(x_min, label_x),  # ensure it's inside visible range
                y=log2_fc,
                s=f"{fc}",
                color="r",
                verticalalignment="center",
                horizontalalignment="left",
                fontsize=45,
            )


    # Axis and legend
    ax.set_xticks(group_centers)
    ax.set_xticklabels(names, fontsize=45)
    ax.set_ylabel(y_label, fontsize=45)
    ax.tick_params(axis="x", labelsize=45)
    ax.tick_params(axis="y", labelsize=45)
    plt.setp(ax.get_xticklabels(), rotation=75, ha="right")

    method_labels = [
        "CruxLFQ", "FlashLFQ", "Sage", "MaxQuant", "IonQuant", "ProteomicsLFQ"
    ]
   
    handles = [
        plt.Line2D([0], [0], color=colors[i], lw=10, label=method_labels[i])
        for i in range(len(colors))
    ]
    ax.legend(handles=handles, loc="best", fontsize=45, frameon=False)

    ax.set_ylim(-3, 2)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)

# --- Human Data ---
human_data_prepped = [
    preprocess(HUMAN_CRUX),
    preprocess(HUMAN_FLASH),
    preprocess(HUMAN_SAGE),
    preprocess(HUMAN_MAX_MAXQUANT),
    preprocess(HUMAN_ION),
    preprocess(HUMAN_PROTEOMICS),
]

plot_combined_boxplot(
    human_data_prepped,
    title="Human",
    filename="human_merged_boxplot.pdf",
)

# --- E. coli Data ---
ecoli_data_prepped = [
    preprocess(ECOLI_CRUX),
    preprocess(ECOLI_FLASH),
    preprocess(ECOLI_SAGE),
    preprocess(ECOLI_MAX_MAXQUANT),
    preprocess(ECOLI_ION),
    preprocess(ECOLI_PROTEOMICS),
]

plot_combined_boxplot(
    ecoli_data_prepped,
    title="E. coli",
    filename="ecoli_merged_boxplot.pdf",
)
