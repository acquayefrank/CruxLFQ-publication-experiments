import numpy as np
import pandas as pd
from scipy.stats import pearsonr


base = "../results/formatted_results/"


def load_data(file_path, sep="\t"):
    """Loads data from the given file path."""
    return pd.read_csv(file_path, sep=sep)


# Helper for printing correlation stats
def print_corr_stats(name, correlations):
    average_correlation = np.mean(correlations)
    std_dev_correlation = np.std(correlations)
    print(name)
    print(f"Pearson correlations: {correlations}")
    print(f"Average correlation: {average_correlation:.5f}")
    print(f"Standard deviation of correlations: {std_dev_correlation:.5f}")


# File paths
def get_files():
    return {
        "crux": base + "crux-lfq-mod-pep.txt_formatted",
        "flash": base + "FlashLFQ+mods+protein_id_modpep.txt_formatted",
        "ion": base + "ionquant_combined_peptide.tsv_formatted",
        "maxq": base + "maxquant_peptides.txt_formatted",
        "sage": base + "sage_lfq.tsv_formatted",
        "proteomics": base + "proteomicslfq.mzTab_formatted",
    }


files = get_files()
df_crux = load_data(files["crux"])
df_flash = load_data(files["flash"])
df_ion = load_data(files["ion"])
df_max = load_data(files["maxq"])
df_sage = load_data(files["sage"])
df_proteomics = load_data(files["proteomics"])

# Column mappings for each tool
crux_cols = [
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_06_161103_A1_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_07_161103_A2_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_16_161103_A3_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_17_161103_A4_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_08_161103_B2_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_15_161103_B3_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_18_161103_B4_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_24_161103_C1_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_09_161103_C2_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_14_161103_C3_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_19_161103_C4_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_22_161103_D1_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_11_161103_D2_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_12_161103_D3_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_21_161103_D4_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_23_161103_E1_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_10_161103_E2_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_13_161103_E3_HCD_OT_4ul.mzML",
    "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_20_161103_E4_HCD_OT_4ul.mzML",
]
flash_cols = [
    "Intensity_B02_06_161103_A1_HCD_OT_4ul",
    "Intensity_B02_07_161103_A2_HCD_OT_4ul",
    "Intensity_B02_16_161103_A3_HCD_OT_4ul",
    "Intensity_B02_17_161103_A4_HCD_OT_4ul",
    "Intensity_B02_001_161103_B1_HCD_OT_4ul",
    "Intensity_B02_08_161103_B2_HCD_OT_4ul",
    "Intensity_B02_15_161103_B3_HCD_OT_4ul",
    "Intensity_B02_18_161103_B4_HCD_OT_4ul",
    "Intensity_B02_24_161103_C1_HCD_OT_4ul",
    "Intensity_B02_09_161103_C2_HCD_OT_4ul",
    "Intensity_B02_14_161103_C3_HCD_OT_4ul",
    "Intensity_B02_19_161103_C4_HCD_OT_4ul",
    "Intensity_B02_22_161103_D1_HCD_OT_4ul",
    "Intensity_B02_11_161103_D2_HCD_OT_4ul",
    "Intensity_B02_12_161103_D3_HCD_OT_4ul",
    "Intensity_B02_21_161103_D4_HCD_OT_4ul",
    "Intensity_B02_23_161103_E1_HCD_OT_4ul",
    "Intensity_B02_10_161103_E2_HCD_OT_4ul",
    "Intensity_B02_13_161103_E3_HCD_OT_4ul",
    "Intensity_B02_20_161103_E4_HCD_OT_4ul",
]
ion_cols = [
    "Exp1_1 Intensity",
    "Exp1_2 Intensity",
    "Exp1_3 Intensity",
    "Exp1_4 Intensity",
    "Exp2_1 Intensity",
    "Exp2_2 Intensity",
    "Exp2_3 Intensity",
    "Exp2_4 Intensity",
    "Exp3_1 Intensity",
    "Exp3_2 Intensity",
    "Exp3_3 Intensity",
    "Exp3_4 Intensity",
    "Exp4_1 Intensity",
    "Exp4_2 Intensity",
    "Exp4_3 Intensity",
    "Exp4_4 Intensity",
    "Exp5_1 Intensity",
    "Exp5_2 Intensity",
    "Exp5_3 Intensity",
    "Exp5_4 Intensity"

]
maxq_cols = [
    "Intensity A1",
    "Intensity A2",
    "Intensity A3",
    "Intensity A4",
    "Intensity B1",
    "Intensity B2",
    "Intensity B3",
    "Intensity B4",
    "Intensity C1",
    "Intensity C2",
    "Intensity C3",
    "Intensity C4",
    "Intensity D1",
    "Intensity D2",
    "Intensity D3",
    "Intensity D4",
    "Intensity E1",
    "Intensity E2",
    "Intensity E3",
    "Intensity E4",
]
sage_cols = [
    "B02_06_161103_A1_HCD_OT_4ul.mzML",
    "B02_07_161103_A2_HCD_OT_4ul.mzML",
    "B02_16_161103_A3_HCD_OT_4ul.mzML",
    "B02_17_161103_A4_HCD_OT_4ul.mzML",
    "B02_001_161103_B1_HCD_OT_4ul.mzML",
    "B02_08_161103_B2_HCD_OT_4ul.mzML",
    "B02_15_161103_B3_HCD_OT_4ul.mzML",
    "B02_18_161103_B4_HCD_OT_4ul.mzML",
    "B02_24_161103_C1_HCD_OT_4ul.mzML",
    "B02_09_161103_C2_HCD_OT_4ul.mzML",
    "B02_14_161103_C3_HCD_OT_4ul.mzML",
    "B02_19_161103_C4_HCD_OT_4ul.mzML",
    "B02_22_161103_D1_HCD_OT_4ul.mzML",
    "B02_11_161103_D2_HCD_OT_4ul.mzML",
    "B02_12_161103_D3_HCD_OT_4ul.mzML",
    "B02_21_161103_D4_HCD_OT_4ul.mzML",
    "B02_23_161103_E1_HCD_OT_4ul.mzML",
    "B02_10_161103_E2_HCD_OT_4ul.mzML",
    "B02_13_161103_E3_HCD_OT_4ul.mzML",
    "B02_20_161103_E4_HCD_OT_4ul.mzML",
]
proteomics_cols = [
    "A1",
    "A2",
    "A3",
    "A4",
    "B1",
    "B2",
    "B3",
    "B4",
    "C1",
    "C2",
    "C3",
    "C4",
    "D1",
    "D2",
    "D3",
    "D4",
    "E1",
    "E2",
    "E3",
    "E4",
]

# 1. Crux vs Flash
print("******** Crux vs Flash ********")
df = pd.merge(df_crux, df_flash, on="id", how="inner")
correlations = []
for col_crux, col_flash in zip(crux_cols, flash_cols):
    df_temp = df[[col_crux, col_flash]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_crux]), np.log2(df_temp[col_flash]))
        correlations.append(corr)
print_corr_stats("Crux vs Flash", correlations)

# 2. Crux vs IonQuant
print("******** Crux vs IonQuant ********")
df = pd.merge(df_crux, df_ion, on="id", how="inner")
correlations = []
for col_crux, col_ion in zip(crux_cols, ion_cols):
    df_temp = df[[col_crux, col_ion]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_crux]), np.log2(df_temp[col_ion]))
        correlations.append(corr)
print_corr_stats("Crux vs IonQuant", correlations)

# 3. Crux vs MaxQuant
print("******** Crux vs MaxQuant ********")
df = pd.merge(df_crux, df_max, on="id", how="inner")
correlations = []
for col_crux, col_max in zip(crux_cols, maxq_cols):
    df_temp = df[[col_crux, col_max]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_crux]), np.log2(df_temp[col_max]))
        correlations.append(corr)
print_corr_stats("Crux vs MaxQuant", correlations)

# 4. Crux vs Sage
print("******** Crux vs Sage ********")
df = pd.merge(df_crux, df_sage, on="id", how="inner")
correlations = []
for col_crux, col_sage in zip(crux_cols, sage_cols):
    df_temp = df[[col_crux, col_sage]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_crux]), np.log2(df_temp[col_sage]))
        correlations.append(corr)
print_corr_stats("Crux vs Sage", correlations)

# 5. Flash vs IonQuant
print("******** Flash vs IonQuant ********")
df = pd.merge(df_flash, df_ion, on="id", how="inner")
correlations = []
for col_flash, col_ion in zip(flash_cols, ion_cols):
    df_temp = df[[col_flash, col_ion]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_flash]), np.log2(df_temp[col_ion]))
        correlations.append(corr)
print_corr_stats("Flash vs IonQuant", correlations)

# 6. Flash vs MaxQuant
print("******** Flash vs MaxQuant ********")
df = pd.merge(df_flash, df_max, on="id", how="inner")
correlations = []
for col_flash, col_max in zip(flash_cols, maxq_cols):
    df_temp = df[[col_flash, col_max]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_flash]), np.log2(df_temp[col_max]))
        correlations.append(corr)
print_corr_stats("Flash vs MaxQuant", correlations)

# 7. Flash vs Sage
print("******** Flash vs Sage ********")
df = pd.merge(df_flash, df_sage, on="id", how="inner")
correlations = []
for col_flash, col_sage in zip(flash_cols, sage_cols):
    df_temp = df[[col_flash, col_sage]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_flash]), np.log2(df_temp[col_sage]))
        correlations.append(corr)
print_corr_stats("Flash vs Sage", correlations)

# 8. IonQuant vs MaxQuant
print("******** IonQuant vs MaxQuant ********")
df = pd.merge(df_ion, df_max, on="id", how="inner")
correlations = []
for col_ion, col_max in zip(ion_cols, maxq_cols):
    df_temp = df[[col_ion, col_max]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_ion]), np.log2(df_temp[col_max]))
        correlations.append(corr)
print_corr_stats("IonQuant vs MaxQuant", correlations)

# 9. IonQuant vs Sage
print("******** IonQuant vs Sage ********")
df = pd.merge(df_ion, df_sage, on="id", how="inner")
correlations = []
for col_ion, col_sage in zip(ion_cols, sage_cols):
    df_temp = df[[col_ion, col_sage]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_ion]), np.log2(df_temp[col_sage]))
        correlations.append(corr)
print_corr_stats("IonQuant vs Sage", correlations)

# 10. MaxQuant vs Sage
print("******** MaxQuant vs Sage ********")
df = pd.merge(df_max, df_sage, on="id", how="inner")
correlations = []
for col_max, col_sage in zip(maxq_cols, sage_cols):
    df_temp = df[[col_max, col_sage]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_max]), np.log2(df_temp[col_sage]))
        correlations.append(corr)
print_corr_stats("MaxQuant vs Sage", correlations)

# 11. ProteomicsLFQ vs MaxQuant
print("******** ProteomicsLFQ vs MaxQuant ********")
df = pd.merge(df_proteomics, df_max, on="id", how="inner")
correlations = []
for col_prot, col_max in zip(proteomics_cols, maxq_cols):
    df_temp = df[[col_prot, col_max]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_prot]), np.log2(df_temp[col_max]))
        correlations.append(corr)
print_corr_stats("ProteomicsLFQ vs MaxQuant", correlations)

# 12. ProteomicsLFQ vs Sage
print("******** ProteomicsLFQ vs Sage ********")
df = pd.merge(df_proteomics, df_sage, on="id", how="inner")
correlations = []
for col_prot, col_sage in zip(proteomics_cols, sage_cols):
    df_temp = df[[col_prot, col_sage]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_prot]), np.log2(df_temp[col_sage]))
        correlations.append(corr)
print_corr_stats("ProteomicsLFQ vs Sage", correlations)

# 13. ProteomicsLFQ vs IonQuant
print("******** ProteomicsLFQ vs IonQuant ********")
df = pd.merge(df_proteomics, df_ion, on="id", how="inner")
correlations = []
for col_prot, col_ion in zip(proteomics_cols, ion_cols):
    df_temp = df[[col_prot, col_ion]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_prot]), np.log2(df_temp[col_ion]))
        correlations.append(corr)
print_corr_stats("ProteomicsLFQ vs IonQuant", correlations)

# 14. ProteomicsLFQ vs Flash
print("******** ProteomicsLFQ vs Flash ********")
df = pd.merge(df_proteomics, df_flash, on="id", how="inner")
correlations = []
for col_prot, col_flash in zip(proteomics_cols, flash_cols):
    df_temp = df[[col_prot, col_flash]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_prot]), np.log2(df_temp[col_flash]))
        correlations.append(corr)
print_corr_stats("ProteomicsLFQ vs Flash", correlations)

# 15. ProteomicsLFQ vs Crux
print("******** ProteomicsLFQ vs Crux ********")
df = pd.merge(df_proteomics, df_crux, on="id", how="inner")
correlations = []
for col_prot, col_crux in zip(proteomics_cols, crux_cols):
    df_temp = df[[col_prot, col_crux]].replace(0, np.nan).dropna()
    if not df_temp.empty:
        corr, _ = pearsonr(np.log2(df_temp[col_prot]), np.log2(df_temp[col_crux]))
        correlations.append(corr)
print_corr_stats("ProteomicsLFQ vs Crux", correlations)
