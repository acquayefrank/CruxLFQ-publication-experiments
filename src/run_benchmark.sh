crux=$(realpath "../bin/crux-lfq")
flash_lfq=$(realpath "../bin/FlashLFQWrapper")

results_folder=$(realpath "../results")
timestamp=$(date +%F_%H%M%S)

flash_output_dir="${results_folder}/${timestamp}_FlashLFQ-benchmark"
crux_output_dir="${results_folder}/${timestamp}_crux-lfq-benchmark"

# The code below assumes that the crux_to_flash results file is named "all_identifications.tsv"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the FlashLFQ command.
identifications=$(find "$results_folder" -type f -name "all_identifications.tsv" -exec realpath {} \; | head -n 1)

spectrum_files_folder=$(realpath "../data/spectrum_files")
extension=".mzML"
# The code below assumes that the percolator target PSMs file is named "percolator.target.psms.txt"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the crux lfq command.
percolator=$(find "$results_folder" -type f -name "percolator.target.psms.txt" -exec realpath {} \; | head -n 1)

for i in {1..3}; do
    shuffled_files=($(printf "%s\n" "${spectra_files[@]}" | shuf))
    echo "Run $i:"
    for num_files in 1 2 4 8 16 20; do
        echo "Selecting $num_files files:"
        selected_files=("${shuffled_files[@]:0:$num_files}")
        
        # Run FlashLFQ
        flash_log_file="${flash_output_dir}/iteration_${i}_${num_files}.log"
        mkdir -p "$flash_output_dir"
        flash_command="$flash_lfq $identifications --out ${flash_output_dir}/FlashLFQ_${num_files}_files --spectrum-files ${selected_files[@]}"
        /usr/bin/time -v $flash_command > "${flash_log_file}" 2>&1
        
        # Run Crux LFQ
        crux_log_file="${crux_output_dir}/iteration_${i}_${num_files}.log"
        mkdir -p "$crux_output_dir"
        crux_command="$crux lfq $percolator ${selected_files[@]} --output-dir ${crux_output_dir}/CruxLFQ_${num_files}_files --spectrum-parser mstoolkit"
        /usr/bin/time -v $crux_command > "${crux_log_file}" 2>&1
    done
done