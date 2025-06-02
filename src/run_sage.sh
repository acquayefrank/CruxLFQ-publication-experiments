#!/bin/bash

# Change the paths below to paths that are reflective of what if on your system.
sage=$(realpath "../bin/sage")
spectrum_files_folder=$(realpath "../data/spectrum_files")
sage_settings=/tear/acquayefrank/sage/sage_settings.json

export SAGE_LOG=trace
export RAYON_NUM_THREADS=4

extension=".mzML"
results_folder=$(realpath "../results")
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_sage-output"


nohup $sage $sage_settings $spectrum_files_folder/*$extension -o $output_dir --batch-size 1 > nohup_sage.txt 2>&1 &
