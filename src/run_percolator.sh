#!/bin/bash


# This is a Bash Script that runs percolator
crux=$(realpath "../bin/crux-percolator")
results_folder=$(realpath "../results")


# The code below assumes that the make-pin results file is named "make-pin.pin"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the percolator command.

make_pin=$(find "$results_folder" -type f -name "make-pin.pin" -exec realpath {} \; | head -n 1)
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_crux-percolator-output"
fixed_params="--output-dir $output_dir"


nohup $crux percolator $make_pin $fixed_params > nohup_percolator.txt &
echo "Percolator has been started. Check nohup_percolator.txt for progress and output in $output_dir."
# Note: Ensure that the paths to crux and make-pin.pin are correct.
# End of script