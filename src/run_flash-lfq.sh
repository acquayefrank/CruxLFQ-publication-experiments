flash_lfq=$(realpath "../bin/FlashLFQWrapperReplicates")
results_folder=$(realpath "../results")
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_FlashLFQ-output"

# The code below assumes that the crux_to_flash results file is named "all_identifications.tsv"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the FlashLFQ command.
identifications=$(find "$results_folder" -type f -name "all_identifications.tsv" -exec realpath {} \; | head -n 1)


nohup $flash_lfq $identifications --out $output_dir > nohup_FlashLFQ.txt &