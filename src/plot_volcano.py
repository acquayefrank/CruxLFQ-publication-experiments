import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# --- Data loading and merging (same as plot_fig2.py) ---
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

# --- Volcano plot function ---
def volcano_data_fast(df, cols_map, groupA_keys, groupB_keys):
    groupA = [cols_map[k] for k in groupA_keys if cols_map[k] in df.columns]
    groupB = [cols_map[k] for k in groupB_keys if cols_map[k] in df.columns]
    dataA = df[groupA].replace(0, np.nan)
    dataB = df[groupB].replace(0, np.nan)
    # Only keep rows with at least one non-NaN in each group
    mask = dataA.notna().all(axis=1) & dataB.notna().all(axis=1)
    dataA = dataA[mask]
    dataB = dataB[mask]
    # log2 transform
    logA = np.log2(dataA)
    logB = np.log2(dataB)
    meanA = logA.mean(axis=1)
    meanB = logB.mean(axis=1)
    log2fc = meanA - meanB
    # Vectorized t-test
    tstat, pvals = ttest_ind(logA.values, logB.values, axis=1, equal_var=False, nan_policy='omit')
    return log2fc, -np.log10(pvals)

methods = [
    ('CruxLFQ', crux_cols_map),
    ('FlashLFQ', flash_cols_map),
    ('MaxQuant', maxquant_cols_map),
    ('IonQuant', ion_cols_map),
    ('ProteomicsLFQ', proteomics_cols_map)
]
colors = [
    'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple'
]
group_pairs = [
    (['A1', 'A2', 'A3', 'A4'], ['B1', 'B2', 'B3', 'B4'], 'A/B'),
    (['A1', 'A2', 'A3', 'A4'], ['C1', 'C2', 'C3', 'C4'], 'A/C'),
    (['A1', 'A2', 'A3', 'A4'], ['D1', 'D2', 'D3', 'D4'], 'A/D'),
    (['A1', 'A2', 'A3', 'A4'], ['E1', 'E2', 'E3', 'E4'], 'A/E'),
]

fig, axes = plt.subplots(len(group_pairs), len(methods), figsize=(6*len(methods), 5*len(group_pairs)), sharey='row')
if len(methods) == 1:
    axes = axes[:, None]
if len(group_pairs) == 1:
    axes = axes[None, :]
for row, (groupA, groupB, label) in enumerate(group_pairs):
    for col, ((method, cmap), color) in enumerate(zip(methods, colors)):
        ax = axes[row, col]
        log2fc, neglog10p = volcano_data_fast(df, cmap, groupA, groupB)
        ax.scatter(log2fc, neglog10p, alpha=0.4, s=18, color=color)
        ax.set_xlabel('log2 fold change (%s)' % label)
        ax.set_title(f'{method}')
        ax.axvline(x=1, color='grey', linestyle='--', lw=1)
        ax.axvline(x=-1, color='grey', linestyle='--', lw=1)
        ax.axhline(y=2, color='grey', linestyle='--', lw=1)
    axes[row,0].set_ylabel('-log10(p-value)\n%s' % label)
plt.tight_layout()
plt.savefig('volcano.pdf')
plt.close()