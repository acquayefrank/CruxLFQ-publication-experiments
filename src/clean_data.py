import pandas as pd
import re
import os
import traceback

def format_sage_lfq(input_file, out_file, q_value_cutoff=0.01):
    df = pd.read_csv(input_file, sep='\t')
    df['id'] = df['peptide'].str.replace(r'\[\+57\.0215\]', '', regex=True)
    df['id'] = df['id'].str.replace(r'M\[\+15\.994915\]', 'm', regex=True)
    df = df[df['q_value_B02_001_161103_B1_HCD_OT_4ul'] <= q_value_cutoff]

    def fix_n_term_mod(x):
        if '[+42.010567]' in x:
            x = x.replace('[+42.010567]-', '')
            chars = list(x)
            chars[0] = 'n-' + chars[0]
            return ''.join(chars)
        return x

    df['id'] = df['id'].apply(fix_n_term_mod)
    df.to_csv(out_file, sep='\t', index=False)


def format_crux_lfq(input_file, out_file, q_value_cutoff=0.01):
    df = pd.read_csv(input_file, sep='\t')
    df['id'] = df['Sequence'].str.replace(r'M\[15\.9949\]', 'm', regex=True)

    def fix_n_term_mod(x):
        if '[42.0106]' in x:
            x = x.replace('[42.0106]', '')
            chars = list(x)
            chars[0] = 'n-' + chars[0]
            return ''.join(chars)
        return x

    df['id'] = df['id'].apply(fix_n_term_mod)
    df.to_csv(out_file, sep='\t', index=False)


def format_flash_lfq(input_file, out_file, q_value_cutoff=0.01):
    df = pd.read_csv(input_file, sep='\t')
    df['id'] = df['Sequence'].str.replace(r'M\[15\.9949\]', 'm', regex=True)

    def fix_n_term_mod(x):
        if '[42.0106]' in x:
            x = x.replace('[42.0106]', '')
            chars = list(x)
            chars[0] = 'n-' + chars[0]
            return ''.join(chars)
        return x

    df['id'] = df['id'].apply(fix_n_term_mod)
    df.to_csv(out_file, sep='\t', index=False)


def format_maxquant_lfq(input_file, out_file, q_value_cutoff=0.01):
    df = pd.read_csv(input_file, sep='\t')
    df['id'] = df['Sequence']
    df.to_csv(out_file, sep='\t', index=False)


def format_ionquant_lfq(input_file, out_file, q_value_cutoff=0.01):
    df = pd.read_csv(input_file, sep='\t')
    df['id'] = df['Peptide Sequence'].str.replace(r'\[57\.0215\]', '', regex=True)
    df['id'] = df['id'].str.replace(r'M\[15\.9949\]', 'm', regex=True)

    def fix_n_term_mod(x):
        if 'n[42.0106]' in x:
            x = x.replace('n[42.0106]', '')
            chars = list(x)
            chars[0] = 'n-' + chars[0]
            return ''.join(chars)
        return x

    df['id'] = df['id'].apply(fix_n_term_mod)
    df.to_csv(out_file, sep='\t', index=False)


def process_all_files(folder_path="../results/organized_results"):
    files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

    # Ensure the output directory exists
    formatted_dir = os.path.join(os.path.dirname(folder_path), "formatted_results")
    os.makedirs(formatted_dir, exist_ok=True)

    for file in files:
        file_name = os.path.basename(file)
        output_file = os.path.join(formatted_dir, f"{file_name}_formatted")

        try:
            if "sage" in file_name.lower():
                format_sage_lfq(file, output_file, q_value_cutoff=1)
            elif "crux" in file_name.lower():
                format_crux_lfq(file, output_file)
            elif "flash" in file_name.lower():
                format_flash_lfq(file, output_file)
            elif "maxquant" in file_name.lower():
                format_maxquant_lfq(file, output_file)
            elif "ionquant" in file_name.lower():
                format_ionquant_lfq(file, output_file)
            else:
                print(f"Unknown file type: {file_name}. Skipping.")
                continue
        except Exception as e:
            print(f"Failed to process: {file}")
            print(f"Error: {str(e)}")
            traceback.print_exc()

if __name__ == "__main__":
    process_all_files()
    print("Data formatting completed.")