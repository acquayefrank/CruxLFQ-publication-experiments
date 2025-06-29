import pandas as pd 
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

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

# Extract sets of IDs from each dataframe
sets = {
    'Crux': set(df_crux['id']),
    'FlashLFQ': set(df_flash['id']),
    'Sage': set(df_sage['id']),
    'IonQuant': set(df_ion['id']),
    'MaxQuant': set(df_max_maxquant['id']),
    'ProteomicsLFQ': set(df_proteomics['id'])
}


# Build a DataFrame for UpSet plot
all_ids = set.union(*sets.values())
id_to_membership = {}
for id_ in all_ids:
    membership = tuple(name for name, s in sets.items() if id_ in s)
    id_to_membership.setdefault(membership, 0)
    id_to_membership[membership] += 1

upset_data = from_memberships(id_to_membership.keys(), data=list(id_to_membership.values()))

# Plot
plt.figure(figsize=(10, 6))
UpSet(upset_data, show_counts=True, sort_by='cardinality').plot()
plt.title('UpSet Plot of Peptides Across Tools')
plt.savefig('upset_plot.pdf', dpi=300, bbox_inches='tight')
plt.close()