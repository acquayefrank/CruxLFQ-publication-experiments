import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde, pearsonr


def plot_density_scatter(
    x, y, xlabel, ylabel, out_file, xlim=(10, 40), ylim=(10, 40), fontsize=30
):
    # Convert to numpy arrays to avoid pandas indexing issues
    x = np.asarray(x)
    y = np.asarray(y)
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    pearson_corr, _ = pearsonr(x, y)
    corr_text = f"Pearson Correlation: {pearson_corr:.5f}"
    num_points = len(x)
    points_text = f"n: {num_points}"

    fig, ax = plt.subplots()
    ax.text(
        0.05, 0.95, corr_text, transform=ax.transAxes, fontsize=fontsize,
        verticalalignment="top", horizontalalignment="left",
    )
    ax.text(
        0.95, 0.05, points_text, transform=ax.transAxes, fontsize=fontsize,
        verticalalignment="bottom", horizontalalignment="right",
    )

    x_line = np.linspace(xlim[0], xlim[1], 100)
    plt.plot(x_line, x_line, color="red", label="y=x")
    scatter = ax.scatter(x, y, c=z, s=1)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    fig.set_size_inches(10, 8)
    ax.tick_params(axis="both", which="major", labelsize=20)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    plt.savefig(out_file)
    plt.close(fig)


def main():
    # File paths
    base = "../results/formatted_results/"
    crux = base + "crux-lfq-mod-pep.txt_formatted"
    flash = base + "FlashLFQ+mods+protein_id_modpep.txt_formatted"
    ion = base + "ionquant_combined_peptide.tsv_formatted"
    maxq = base + "maxquant_peptides.txt_formatted"
    sage = base + "sage_lfq.tsv_formatted"
    proteomics = base + "proteomicslfq.mzTab_formatted"

    # 1. Crux vs Flash
    df_crux = pd.read_csv(crux, sep="\t")
    df_flash = pd.read_csv(flash, sep="\t")
    df1 = pd.merge(
        df_crux.sort_values("id"),
        df_flash.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML"
    col_y = "Intensity_B02_001_161103_B1_HCD_OT_4ul"
    df1[col_x] = df1[col_x].replace(0, np.nan)
    df1[col_y] = df1[col_y].replace(0, np.nan)
    df1 = df1.dropna(subset=[col_x, col_y])
    x = np.log2(df1[col_x])
    y = np.log2(df1[col_y])
    plot_density_scatter(
        x, y,
        "CruxLFQ log2 Intensity", "FlashLFQ log2 Intensity",
        "CruxLFQ_vs_FlashLFQ_density_scatter_plot.pdf"
    )

    # 2. Crux vs IonQuant
    df_ion = pd.read_csv(ion, sep="\t")
    df2 = pd.merge(
        df_crux.sort_values("id"),
        df_ion.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML"
    col_y = "Exp2_1 Intensity"
    df2[col_x] = df2[col_x].replace(0, np.nan)
    df2[col_y] = df2[col_y].replace(0, np.nan)
    df2 = df2.dropna(subset=[col_x, col_y])
    x = np.log2(df2[col_x])
    y = np.log2(df2[col_y])
    plot_density_scatter(
        x, y,
        "CruxLFQ log2 Intensity", "IonQuant log2 Intensity",
        "CruxLFQ_vs_IonQuant_density_scatter_plot.pdf"
    )

    # 3. Crux vs MaxQuant
    df_max = pd.read_csv(maxq, sep="\t")
    df3 = pd.merge(
        df_crux.sort_values("id"),
        df_max.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML"
    col_y = "Intensity B1"
    df3[col_x] = df3[col_x].replace(0, np.nan)
    df3[col_y] = df3[col_y].replace(0, np.nan)
    df3 = df3.dropna(subset=[col_x, col_y])
    x = np.log2(df3[col_x])
    y = np.log2(df3[col_y])
    plot_density_scatter(
        x, y,
        "CruxLFQ log2 Intensity", "MaxQuant log2 Intensity",
        "CruxLFQ_vs_MaxLFQ_density_scatter_plot.pdf"
    )

    # 4. Crux vs Sage
    df_sage = pd.read_csv(sage, sep="\t")
    df4 = pd.merge(
        df_crux.sort_values("id"),
        df_sage.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML"
    col_y = "B02_001_161103_B1_HCD_OT_4ul.mzML"
    df4[col_x] = df4[col_x].replace(0, np.nan)
    df4[col_y] = df4[col_y].replace(0, np.nan)
    df4 = df4.dropna(subset=[col_x, col_y])
    x = np.log2(df4[col_x])
    y = np.log2(df4[col_y])
    plot_density_scatter(
        x, y,
        "CruxLFQ log2 Intensity", "Sage log2 Intensity",
        "CruxLFQ_vs_SageLFQ_density_scatter_plot.pdf"
    )

    # 5. Flash vs IonQuant
    df5 = pd.merge(
        df_flash.sort_values("id"),
        df_ion.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_B02_001_161103_B1_HCD_OT_4ul"
    col_y = "Exp2_1 Intensity"
    df5[col_x] = df5[col_x].replace(0, np.nan)
    df5[col_y] = df5[col_y].replace(0, np.nan)
    df5 = df5.dropna(subset=[col_x, col_y])
    x = np.log2(df5[col_x])
    y = np.log2(df5[col_y])
    plot_density_scatter(
        x, y,
        "FlashLFQ log2 Intensity", "IonQuant log2 Intensity",
        "FlashLFQ_vs_IonQuant_density_scatter_plot.pdf"
    )

    # 6. Flash vs MaxQuant
    df6 = pd.merge(
        df_flash.sort_values("id"),
        df_max.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_B02_001_161103_B1_HCD_OT_4ul"
    col_y = "Intensity B1"
    df6[col_x] = df6[col_x].replace(0, np.nan)
    df6[col_y] = df6[col_y].replace(0, np.nan)
    df6 = df6.dropna(subset=[col_x, col_y])
    x = np.log2(df6[col_x])
    y = np.log2(df6[col_y])
    plot_density_scatter(
        x, y,
        "FlashLFQ log2 Intensity", "MaxQuant log2 Intensity",
        "FlashLFQ_vs_MaxLFQ_density_scatter_plot.pdf"
    )

    # 7. Flash vs Sage
    df7 = pd.merge(
        df_flash.sort_values("id"),
        df_sage.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity_B02_001_161103_B1_HCD_OT_4ul"
    col_y = "B02_001_161103_B1_HCD_OT_4ul.mzML"
    df7[col_x] = df7[col_x].replace(0, np.nan)
    df7[col_y] = df7[col_y].replace(0, np.nan)
    df7 = df7.dropna(subset=[col_x, col_y])
    x = np.log2(df7[col_x])
    y = np.log2(df7[col_y])
    plot_density_scatter(
        x, y,
        "FlashLFQ log2 Intensity", "Sage log2 Intensity",
        "FlashLFQ_vs_Sage_density_scatter_plot.pdf"
    )

    # 8. IonQuant vs MaxQuant
    df8 = pd.merge(
        df_ion.sort_values("id"),
        df_max.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Exp2_1 Intensity"
    col_y = "Intensity B1"
    df8[col_x] = df8[col_x].replace(0, np.nan)
    df8[col_y] = df8[col_y].replace(0, np.nan)
    df8 = df8.dropna(subset=[col_x, col_y])
    x = np.log2(df8[col_x])
    y = np.log2(df8[col_y])
    plot_density_scatter(
        x, y,
        "IonQuant log2 Intensity", "MaxQuant log2 Intensity",
        "MaxLFQ_vs_IonQuant_density_scatter_plot.pdf"
    )

    # 9. IonQuant vs Sage
    df9 = pd.merge(
        df_ion.sort_values("id"),
        df_sage.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Exp2_1 Intensity"
    col_y = "B02_001_161103_B1_HCD_OT_4ul.mzML"
    df9[col_x] = df9[col_x].replace(0, np.nan)
    df9[col_y] = df9[col_y].replace(0, np.nan)
    df9 = df9.dropna(subset=[col_x, col_y])
    x = np.log2(df9[col_x])
    y = np.log2(df9[col_y])
    plot_density_scatter(
        y, x,
        "Sage log2 Intensity", "IonQuant log2 Intensity",
        "Sage_vs_IonQuant_density_scatter_plot.pdf"
    )

    # 10. MaxQuant vs Sage
    df10 = pd.merge(
        df_max.sort_values("id"),
        df_sage.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "Intensity B1"
    col_y = "B02_001_161103_B1_HCD_OT_4ul.mzML"
    df10[col_x] = df10[col_x].replace(0, np.nan)
    df10[col_y] = df10[col_y].replace(0, np.nan)
    df10 = df10.dropna(subset=[col_x, col_y])
    x = np.log2(df10[col_x])
    y = np.log2(df10[col_y])
    plot_density_scatter(
        y, x,
        "Sage log2 Intensity", "MaxQuant log2 Intensity", 
        "MaxLFQ_vs_Sage_density_scatter_plot.pdf"
    )

    # 11. ProteomicsLFQ vs MaxQuant
    df_proteomics = pd.read_csv(proteomics, sep="\t")
    df11 = pd.merge(
        df_proteomics.sort_values("id"),
        df_max.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "B1"
    col_y = "Intensity B1"
    df11[col_x] = df11[col_x].replace(0, np.nan)
    df11[col_y] = df11[col_y].replace(0, np.nan)
    df11 = df11.dropna(subset=[col_x, col_y])
    x = np.log2(df11[col_x])
    y = np.log2(df11[col_y])
    plot_density_scatter(
        y, x,
        "ProteomicsLFQ log2 Intensity", "MaxQuant log2 Intensity",
        "ProteomicsLFQ_vs_MaxLFQ_density_scatter_plot.pdf"
    )

    # 12. ProteomicsLFQ vs Sage
    df12 = pd.merge(
        df_proteomics.sort_values("id"),
        df_sage.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "B1"
    col_y = "B02_001_161103_B1_HCD_OT_4ul.mzML"
    df12[col_x] = df12[col_x].replace(0, np.nan)
    df12[col_y] = df12[col_y].replace(0, np.nan)
    df12 = df12.dropna(subset=[col_x, col_y])
    x = np.log2(df12[col_x])
    y = np.log2(df12[col_y])
    plot_density_scatter(
        y, x,
        "ProteomicsLFQ log2 Intensity", "Sage log2 Intensity",
        "ProteomicsLFQ_vs_Sage_density_scatter_plot.pdf"
    )

    # 13. ProteomicsLFQ vs IonQuant
    df13 = pd.merge(
        df_proteomics.sort_values("id"),
        df_ion.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "B1"
    col_y = "Exp2_1 Intensity"  
    df13[col_x] = df13[col_x].replace(0, np.nan)
    df13[col_y] = df13[col_y].replace(0, np.nan)
    df13 = df13.dropna(subset=[col_x, col_y])
    x = np.log2(df13[col_x])
    y = np.log2(df13[col_y])
    plot_density_scatter(
        y, x,
        "IonQuant log2 Intensity", "ProteomicsLFQ log2 Intensity",
        "IonQuant_vs_ProteomicsLFQ_density_scatter_plot.pdf"
    )

    # 14. ProteomicsLFQ vs Flash
    df14 = pd.merge(
        df_proteomics.sort_values("id"),
        df_flash.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "B1"
    col_y = "Intensity_B02_001_161103_B1_HCD_OT_4ul"
    df14[col_x] = df14[col_x].replace(0, np.nan)
    df14[col_y] = df14[col_y].replace(0, np.nan)
    df14 = df14.dropna(subset=[col_x, col_y])
    x = np.log2(df14[col_x])
    y = np.log2(df14[col_y])
    plot_density_scatter(
        y, x,
        "Flash log2 Intensity", "ProteomicsLFQ log2 Intensity",
        "Flash_vs_ProteomicsLFQ_density_scatter_plot.pdf"
    )

    # 15. ProteomicsLFQ vs Crux
    df15 = pd.merge(
        df_proteomics.sort_values("id"),
        df_crux.sort_values("id"),
        left_on="id", right_on="id", how="inner"
    )
    col_x = "B1"
    col_y = "Intensity_/home/acquayefrank/CruxLFQ-publication-experiments/data/spectrum_files/B02_001_161103_B1_HCD_OT_4ul.mzML"
    df15[col_x] = df15[col_x].replace(0, np.nan)
    df15[col_y] = df15[col_y].replace(0, np.nan)
    df15 = df15.dropna(subset=[col_x, col_y])
    x = np.log2(df15[col_x])
    y = np.log2(df15[col_y])
    plot_density_scatter(
        y, x,
        "Crux log2 Intensity", "ProteomicsLFQ log2 Intensity",
        "Crux_vs_ProteomicsLFQ_density_scatter_plot.pdf"
    )


if __name__ == "__main__":
    main()