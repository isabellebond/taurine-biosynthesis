import pandas as pd
import os
import tarfile
import numpy as np

def extract_coloc_results(folder, outfile):

    all_data = []
    failed_files = []
    # Iterate over all .tar files in the folder
    for filename in os.listdir(folder):
        if filename.endswith('.tar.gz'):
            print(f'extracting {filename}...')
            tar_path = os.path.join(folder, filename)
            try:
                with tarfile.open(tar_path, 'r') as tar:
                    # Look for coloc_results.tsv inside the tar
                    for member in tar.getmembers():
                        if os.path.basename(member.name) == 'coloc_results.tsv.gz':
                            f = tar.extractfile(member)
                            if f is not None:
                                df = pd.read_csv(f, sep='\t', compression = 'gzip')
                                df['source_tarfile'] = filename
                                df['gene'] = filename.split('.')[0]
                                all_data.append(df)
                            break  # Assume only one coloc_results.tsv per tar
            except tarfile.ReadError:
                failed_files.append(filename)

    # Combine all data and save to a single TSV file
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_df.to_csv(outfile, sep='\t', index=False)
        print(f"Combined file written to {outfile}")
    else:
        print("No coloc_results.tsv files found.")
    
    if failed_files:
        with open(f"{outfile}.failed", 'w') as f:
            for item in failed_files:
                f.write(item.strip() + "\n")

    
    return combined_df