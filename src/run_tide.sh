#!/bin/bash


crux=$(realpath "../bin/crux-4.3.Linux.x86_64/bin/crux")
fasta_file=$(realpath "../data/Human_ecoli_trypsin_1501v_uniprot_sprot.fasta")
spectrum_files_folder=$(realpath "../data/spectrum_files")
results_folder=$(realpath "../results")

extension=".mzML"
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_crux_output"
fixed_params="--output-dir $output_dir --concat T --num-threads 6"
index_dir=${results_folder}/HUMAN_ECOLI


echo CRUX: $crux
echo FASTA: $fasta_file
echo SPECTRA: $spectrum_files_folder
echo RESULTS: $results_folder

echo "Starting Tide index and search..."

$crux tide-index --mods-spec 1M+15.9949 --nterm-peptide-mods-spec 1X+42.010565 --missed-cleavages 2 --overwrite T --output-dir $index_dir $fasta_file $index_dir
$crux tide-search --precursor-window 20 --precursor-window-type ppm $fixed_params $spectrum_files_folder/*$extension $index_dir
echo "Tide search completed. Results are in $index_dir and $output_dir"
# End of script
# Note: Ensure that the paths to crux, fasta file, and spectrum files are correct.