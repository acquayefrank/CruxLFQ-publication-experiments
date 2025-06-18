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
df_sage = pd.read_csv(sage, sep="\t")
df_ion = pd.read_csv(ion, sep="\t")
df_max_maxquant = pd.read_csv(maxq, sep="\t")
df_proteomics = pd.read_csv(proteomics, sep="\t")

filtered_df = df_sage[df_sage["q_value_B02_001_161103_B1_HCD_OT_4ul"] <= 1]

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
    filtered_df,
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

def process_data(categories, key, row):
    crux_to_max = {
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
    crux_to_ion = {
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
    crux_to_max_maxquant = {
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

    crux_to_proteomics = {
        "A1": "A1", "A2": "A2", "A3": "A3", "A4": "A4",
        "B1": "B1", "B2": "B2", "B3": "B3", "B4": "B4",
        "C1": "C1", "C2": "C2", "C3": "C3", "C4": "C4",
        "D1": "D1", "D2": "D2", "D3": "D3", "D4": "D4",
        "E1": "E1", "E2": "E2", "E3": "E3", "E4": "E4",
    }


    for category in categories.keys():
        if (
            key == crux_to_max[category]
            or key == crux_to_ion[category]
            or f"_{category}_" in key
            or key == crux_to_max_maxquant[category]
            or key == crux_to_proteomics[category]
        ):
            categories[category].append(row[key])
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

# --- Split data for human and ecoli ---
human_crux_categories = initialize_categories()
human_flash_categories = initialize_categories()
human_sage_categories = initialize_categories()
human_max_categories = initialize_categories()
human_ion_categories = initialize_categories()
human_max_maxquant = initialize_categories()
human_proteomics_categories = initialize_categories()

ecoli_crux_categories = initialize_categories()
ecoli_flash_categories = initialize_categories()
ecoli_sage_categories = initialize_categories()
ecoli_max_categories = initialize_categories()
ecoli_ion_categories = initialize_categories()
ecoli_max_maxquant = initialize_categories()
ecoli_proteomics_categories = initialize_categories()

for index, row in df.iterrows():
    if "HUMAN" in row["Protein ID_df2"]:
        for key in row.index.tolist():
            if "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/" in key:
                human_crux_categories = process_data(human_crux_categories, key, row)
            elif "Intensity_" in key and key.split("_")[-1] == "4ul":
                human_flash_categories = process_data(human_flash_categories, key, row)
            elif ".mzML" in key and '/' not in key:
                human_sage_categories = process_data(
                    human_sage_categories, key, row
                )
            elif "MaxLFQ Intensity" in key and "Exp" in key:
                human_max_categories = process_data(human_max_categories, key, row)
            elif " Intensity" in key:
                human_ion_categories = process_data(human_ion_categories, key, row)
            elif "Intensity " in key:
                human_max_maxquant = process_data(human_max_maxquant, key, row)
            elif key in [
                "A1", "A2", "A3", "A4",
                "B1", "B2", "B3", "B4",
                "C1", "C2", "C3", "C4",
                "D1", "D2", "D3", "D4",
                "E1", "E2", "E3", "E4",
            ]:
                # This is a special case for the categories
                human_proteomics_categories = process_data(human_proteomics_categories, key, row)
    if "ECOLI" in row["Protein ID_df2"]:
        for key in row.index.tolist():
            if "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/" in key:
                ecoli_crux_categories = process_data(ecoli_crux_categories, key, row)
            elif "Intensity_" in key and key.split("_")[-1] == "4ul":
                ecoli_flash_categories = process_data(ecoli_flash_categories, key, row)
            elif ".mzML" in key and '/' not in key:
                ecoli_sage_categories = process_data(
                    ecoli_sage_categories, key, row
                )
            elif "MaxLFQ Intensity" in key and "Exp" in key:
                ecoli_max_categories = process_data(ecoli_max_categories, key, row)
            elif " Intensity" in key:
                ecoli_ion_categories = process_data(ecoli_ion_categories, key, row)
            elif "Intensity " in key:
                ecoli_max_maxquant = process_data(ecoli_max_maxquant, key, row)
            elif key in [
                "A1", "A2", "A3", "A4",
                "B1", "B2", "B3", "B4",
                "C1", "C2", "C3", "C4",
                "D1", "D2", "D3", "D4",
                "E1", "E2", "E3", "E4",
            ]:
                # This is a special case for the categories
                ecoli_proteomics_categories = process_data(ecoli_proteomics_categories, key, row)

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

# --- Plot for HUMAN ---
fig, ax = plt.subplots(figsize=(25, 17))
human_crux_box_data = filter_nan([np.log2(l) for l in HUMAN_CRUX])
human_flash_box_data = filter_nan([np.log2(l) for l in HUMAN_FLASH])
human_sage_box_data = filter_nan([np.log2(l) for l in HUMAN_SAGE])
human_max_box_data = filter_nan([np.log2(l) for l in HUMAN_MAX])
human_ion_box_data = filter_nan([np.log2(l) for l in HUMAN_ION])
human_max_maxquant_box_data = filter_nan([np.log2(l) for l in HUMAN_MAX_MAXQUANT])
human_proteomics_box_data = filter_nan([np.log2(l) for l in HUMAN_PROTEOMICS])


names = ["1.5-fold (B)", "2.0-fold (C)", "2.5-fold (D)", "3.0-fold (E)"]
fold_changes = [1.5, 2.0, 2.5, 3.0]
custom_colors = ["#FFFFFF", "#D3D3D3", "#FFFFFF", "#D3D3D3"]

method_labels = [
    "CruxLFQ", "FlashLFQ", "Sage", "MaxQuant", "IonQuant", "ProteomicsLFQ"
]
colors = [
    "lightblue", "lightgreen", "magenta", "lightyellow", "pink", "lightcoral"
]

# Filter out NaNs and log2 transform
def preprocess(data_tuple):
    return [np.array([np.log2(x) for x in l if not np.isnan(x)]) for l in data_tuple]

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
    for i, color in enumerate(custom_colors):
        ax.axvspan(
            group_centers[i] - spacing * 3,
            group_centers[i] + spacing * 3,
            facecolor=color,
            alpha=0.3,
        )

    # Plot each method
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
        label_x = group_centers[0] - spacing * 5  # inside the first group

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