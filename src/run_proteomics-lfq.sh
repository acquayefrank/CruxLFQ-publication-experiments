#!/bin/bash
set -e

# === Initialize conda in a generic way ===
if ! command -v conda &> /dev/null; then
    echo "❌ conda command not found. Please install Miniconda/Anaconda and ensure 'conda' is in your PATH."
    exit 1
fi

# Initialize conda shell support (for 'conda activate')
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate openms || { echo "❌ Failed to activate conda environment 'openms'"; exit 1; }

echo "✅ Conda environment activated: $(conda info --envs | grep '*' )"

# === CHANGE THESE PATHS ===
DATA_DIR=$(realpath "../data/spectrum_files")
FASTA=$(realpath "../data/Human_ecoli_trypsin_1501v_uniprot_sprot.fasta")
comet=$(realpath "../bin/comet.linux.exe")
base_data_dir=$(realpath "../data")
DESIGN_FILE=$(realpath "./design.tsv")

OUTPUT_DIR=$(realpath "../results/proteomics-lfq")
mkdir -p "$OUTPUT_DIR"

cd "$DATA_DIR"

echo "🔍 Searching for mzML files in $DATA_DIR"
samples=( *.mzML )
if [ ${#samples[@]} -eq 0 ]; then
    echo "❌ No .mzML files found in $DATA_DIR"
    exit 1
fi

# === Step 1: Generate concatenated target-decoy FASTA if not present ===
DECOY_FASTA=${base_data_dir}/Human_ecoli_trypsin_1501v_uniprot_sprot_td.fasta
if [ ! -f "$DECOY_FASTA" ]; then
    echo "🧬 Generating target-decoy FASTA..."
    DecoyDatabase -in "$FASTA" -out "$DECOY_FASTA" -decoy_string DECOY_ -decoy_string_position prefix
else
    echo "🧬 Target-decoy FASTA already exists: $DECOY_FASTA"
fi

# Use the decoy FASTA for all downstream steps
FASTA=$(realpath "$DECOY_FASTA")

# === Step 2: Peptide identification with CometAdapter ===
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.idXML already exists, skipping CometAdapter."
    else
        echo "🔬 Identifying peptides in $sample..."
        CometAdapter \
        -in "$sample" \
        -out "$OUTPUT_DIR/$base.idXML" \
        -database "$FASTA" \
        -comet_executable $comet \
        -missed_cleavages 2 \
        -precursor_mass_tolerance 20.0 \
        -precursor_error_units ppm \
        -variable_modifications "Oxidation (M)" \
        -variable_modifications "Acetyl (N-term)" \
        -enzyme Trypsin \
        -max_variable_mods_in_peptide 3
    fi
done

# === Step 3: Annotate identifications with protein sequence info ===
# This step feels redundant if we are going to run PeptideIndexer later .
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.indexed.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.indexed.idXML already exists, skipping PeptideIndexer."
        elif [ -f "$OUTPUT_DIR/$base.idXML" ]; then
        echo "🔗 Running PeptideIndexer for $base.idXML..."
        PeptideIndexer \
        -in "$OUTPUT_DIR/$base.idXML" \
        -fasta "$FASTA" \
        -out "$OUTPUT_DIR/$base.indexed.idXML"
    else
        echo "⚠️  $OUTPUT_DIR/$base.idXML not found, skipping PeptideIndexer."
    fi
done

# === Step 4: Calculate Posterior Error Probability (PEP) ===
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.pep.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.pep.idXML already exists, skipping IDPosteriorErrorProbability."
        elif [ -f "$OUTPUT_DIR/$base.indexed.idXML" ]; then
        echo "🧮 Calculating Posterior Error Probability for $base.indexed.idXML..."
        IDPosteriorErrorProbability \
        -in "$OUTPUT_DIR/$base.indexed.idXML" \
        -out "$OUTPUT_DIR/$base.pep.idXML"
    else
        echo "⚠️  $OUTPUT_DIR/$base.indexed.idXML not found, skipping PEP calculation."
    fi
done

# === Step 5: Peptide FDR estimation (after PEP) ===
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.fdr.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.fdr.idXML already exists, skipping FDR estimation."
        elif [ -f "$OUTPUT_DIR/$base.pep.idXML" ]; then
        echo "🧮 Estimating peptide FDR for $base.pep.idXML..."
        FalseDiscoveryRate \
        -in "$OUTPUT_DIR/$base.pep.idXML" \
        -out "$OUTPUT_DIR/$base.fdr.idXML" \
        -peptide true \
        -protein false \
        -force
    else
        echo "⚠️  $OUTPUT_DIR/$base.pep.idXML not found, skipping FDR estimation."
    fi
done

# === Step 6: (Optional) Peptide-level filtering (no protein filtering) ===
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.filtered.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.filtered.idXML already exists, skipping IDFilter."
        elif [ -f "$OUTPUT_DIR/$base.fdr.idXML" ]; then
        echo "🧹 Filtering peptides for $base.fdr.idXML... (no protein filtering)"
        IDFilter \
        -in "$OUTPUT_DIR/$base.fdr.idXML" \
        -out "$OUTPUT_DIR/$base.filtered.idXML" \
        -score:peptide 0.01 \
        -threads 4
    else
        echo "⚠️ $OUTPUT_DIR/$base.fdr.idXML not found, skipping IDFilter."
    fi
done

# === Step 6b: Re-annotate filtered IDs with protein sequence info ===
for sample in "${samples[@]}"; do
    base="${sample%.mzML}"
    if [ -f "$OUTPUT_DIR/$base.filtered.indexed.idXML" ]; then
        echo "✅ $OUTPUT_DIR/$base.filtered.indexed.idXML already exists, skipping PeptideIndexer (filtered)."
        elif [ -f "$OUTPUT_DIR/$base.filtered.idXML" ]; then
        echo "🔗 Running PeptideIndexer for $base.filtered.idXML..."
        PeptideIndexer \
        -in "$OUTPUT_DIR/$base.filtered.idXML" \
        -fasta "$FASTA" \
        -out "$OUTPUT_DIR/$base.filtered.indexed.idXML" \
        -write_protein_sequence
    else
        echo "⚠️  $OUTPUT_DIR/$base.filtered.idXML not found, skipping PeptideIndexer (filtered)."
    fi
done

# === Step 8: Peptide quantification with ProteomicsLFQ (per sample) ===
echo "📊 Running peptide-level quantification for each sample..."
ProteomicsLFQ \
-in $(for sample in "${samples[@]}"; do echo "$DATA_DIR/$sample"; done) \
-ids $(for sample in "${samples[@]}"; do echo "$OUTPUT_DIR/${sample%.mzML}.filtered.indexed.idXML"; done) \
-out "$OUTPUT_DIR/lfq.mzTab" \
-out_cxml "$OUTPUT_DIR/lfq.consensusXML" \
-design "$DESIGN_FILE"

echo "✅ Pipeline finished successfully!"
echo "Results saved in: $OUTPUT_DIR"
echo " - lfq.consensusXML"
echo " - lfq.mzTab"
