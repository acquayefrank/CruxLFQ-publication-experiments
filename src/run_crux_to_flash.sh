crux=$(realpath "../bin/crux-to-flash")
spectrum_files_folder=$(realpath "../data/spectrum_files")
results_folder=$(realpath "../results")

extension=".mzML"
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_crux_to_flash-output"

# The code below assumes that the percolator target PSMs file is named "percolator.target.psms.txt"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the crux lfq command.
percolator=$(find "$results_folder" -type f -name "percolator.target.psms.txt" -exec realpath {} \; | head -n 1)

nohup $crux lfq $percolator $spectrum_files_folder/*$extension --output-dir $output_dir --spectrum-parser mstoolkit --mods-spec 1M+15.9949 --nterm-peptide-mods-spec 1X+42.010565  --psm-file-format percolator > nohup_crux_to_flash.txt &
echo "Crux LFQ search started. Results will be saved in $output_dir and nohup output in nohup_crux_to_flash.txt"