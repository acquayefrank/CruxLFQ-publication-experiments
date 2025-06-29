import os
import glob
import pandas as pd

def find_lfq_files(base_dir):
    """
    Find all 'lfq.tsv' files in subfolders containing the word 'sage'.
    
    Args:
        base_dir (str): The base directory to search within.
    
    Returns:
        list: A list of paths to the found 'lfq.tsv' files.
    """
    lfq_files = []
    # Search for directories containing 'sage' in their name
    for root, dirs, files in os.walk(base_dir):
        if 'sage' in root.lower():  # Case-insensitive match for 'sage'
            # Look for 'lfq.tsv' files in the current directory
            lfq_files.extend(glob.glob(os.path.join(root, 'lfq.tsv')))
    return lfq_files

def concat_lfq_files(lfq_files):
    """
    Concatenate multiple 'lfq.tsv' files row-wise, like the R script.
    
    Args:
        lfq_files (list): List of file paths to 'lfq.tsv' files.
    
    Returns:
        pd.DataFrame: A concatenated DataFrame.
    """
    dfs = [pd.read_csv(file, sep='\t') for file in lfq_files]
    return pd.concat(dfs, ignore_index=True)

# Example usage
if __name__ == "__main__":
    base_directory = "../results/"  # Replace with your base directory
    lfq_files = find_lfq_files(base_directory)
    print("Found lfq.tsv files:")
    for file in lfq_files:
        print(file)
    
    if lfq_files:
        concatenated_data = concat_lfq_files(lfq_files)
        print("Concatenated Data:")
        print(concatenated_data.head())  # Display the first few rows of the concatenated DataFrame
        # Optionally, save the concatenated data to a file
        concatenated_data.to_csv("../results/organized_results/sage_lfq.tsv", sep='\t', index=False)