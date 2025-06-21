#!/bin/bash


spectrum_files_folder=$(realpath "../data/spectrum_files")

mapfile -t spectra_files < <(find "$spectra_files_path" -type f -name "*.mzML" -printf '%p\n' | sort -u)

results_folder=$(realpath "../results")
LOG_DIR="${results_folder}/crux_runs_timing_memory_logs"
mkdir -p "$LOG_DIR"

crux=$(realpath "../bin/crux-lfq")
assign_confidence_file=$(realpath "../data/assign-confidence.target.txt")

for i in {1..3}; do
    shuffled_files=($(printf "%s\n" "${spectra_files[@]}" | shuf))
    echo "Run $i:"
    for num_files in 1 2 4 8 16 20; do
        echo "Selecting $num_files files:"
        selected_files=("${shuffled_files[@]:0:$num_files}")
        LOG_FILE="$LOG_DIR/iteration_${i}_${num_files}.log"
        crux="${crux} lfq ${assign_confidence_file} ${selected_files[@]} --overwrite T --spectrum-parser mstoolkit"
        /usr/bin/time -v $crux > "${LOG_FILE}" 2>&1
    done
done