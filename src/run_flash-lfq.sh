flash_lfq=$(realpath "../bin/FlashLFQWrapper")
results_folder=$(realpath "../results")
timestamp=$(date +%F_%H%M%S)
output_dir="${results_folder}/${timestamp}_FlashLFQ-output"
replicate_file=$(realpath "../src/replicate.tsv")

# The code below assumes that the crux_to_flash results file is named "crux-lfq-identifications.txt"
# and is located in the results folder. It finds the first occurrence of this file
# and uses its absolute path for the FlashLFQ command.
# Note if crux-lfq-identifications.txt is not found, try all_identifications.tsv
identifications=$(find "$results_folder" -type f -name "crux-lfq-identifications.txt" -exec realpath {} \; | head -n 1)
if [ -z "$identifications" ]; then
    identifications=$(find "$results_folder" -type f -name "all_identifications.tsv" -exec realpath {} \; | head -n 1)
fi

nohup $flash_lfq $identifications --out $output_dir --replicates $replicate_file > nohup_FlashLFQ.txt &