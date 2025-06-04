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

def merge_lfq_files(lfq_files, key='peptide'):
    """
    Merge multiple 'lfq.tsv' files on a specified key, using parent folder names as suffixes for overlapping columns.
    
    Args:
        lfq_files (list): List of file paths to 'lfq.tsv' files.
        key (str): The column name to merge on.
    
    Returns:
        pd.DataFrame: A merged DataFrame.
    """
    merged_df = None
    for file in lfq_files:
        df = pd.read_csv(file, sep='\t')  # Assuming tab-separated values
        # Extract the parent folder name to use as a suffix
        file_suffix = os.path.basename(os.path.dirname(file))
        print(f"Processing file: {file} with suffix: {file_suffix}")
        # Rename columns to include the suffix (except the key column and columns already containing the suffix)
        df = df.rename(columns={
            col: f"{col}_{file_suffix}" 
            for col in df.columns 
            if col != key and file_suffix not in col
        })
        
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, left_on=f"{key}", right_on=f"{key}", how='outer')
    return merged_df

# Example usage
if __name__ == "__main__":
    base_directory = "../results/"  # Replace with your base directory
    lfq_files = find_lfq_files(base_directory)
    print("Found lfq.tsv files:")
    for file in lfq_files:
        print(file)
    
    if lfq_files:
        merged_data = merge_lfq_files(lfq_files, key='peptide')
        print("Merged Data:")
        print(merged_data.head())  # Display the first few rows of the merged DataFrame
        # Optionally, save the merged data to a file
        merged_data.to_csv("../results/organized_results/sage_lfq.tsv", sep='\t', index=False)