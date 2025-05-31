#!/bin/bash


crux=$(realpath "../bin/crux-lfq")
results_folder=$(realpath "../results")

# The code below assumes that the tide search results file is named "tide-search.txt"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the make-pin command.
tide_search_results=$(find "$results_folder" -type f -name "tide-search.txt" -exec realpath {} \; | head -n 1)

timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_crux-make-pin-output"

fixed_params="--output-dir $output_dir"

nohup $crux make-pin $fixed_params $tide_search_results > nohup_make-pin.txt &
echo "Make-pin has been started. Check the output in $output_dir."
# Note: Ensure that the paths to crux and tide-search.txt are correct.