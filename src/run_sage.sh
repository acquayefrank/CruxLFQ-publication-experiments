#!/bin/bash

sage=$(realpath "../bin/sage")
spectrum_files_folder=$(realpath "../data/spectrum_files")
sage_settings=/tear/acquayefrank/sage/sage_settings.json

export SAGE_LOG=trace
export RAYON_NUM_THREADS=4

extension=".mzML"
results_folder=$(realpath "../results")
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_sage-output"

mkdir -p "$output_dir"

for file in "$spectrum_files_folder"/*$extension; do
    filename=$(basename "$file")
    name="${filename%$extension}"
    file_output_dir="$output_dir/$name"
    mkdir -p "$file_output_dir"
    $sage $sage_settings "$file" -o "$file_output_dir" --batch-size 1
done