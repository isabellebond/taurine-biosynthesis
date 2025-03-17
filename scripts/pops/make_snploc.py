import pandas as pd
import os
import gzip
import shutil
import tempfile
import argparse
import numpy as np

def get_args():
    """
    Parse command-line arguments for extracting information from raw files for GWAS normalization.

    Returns
    -------
    tuple
        A tuple containing:
        - infile (str): The file name of the input file, including the ".gz" extension.
        - outfile (str): The file name of the output file, including the ".gz" extension.

    Notes
    -----
    This function uses argparse to handle command-line arguments. It expects:
    - `--infile`: The path to the input file (must be a gzipped file).
    - `--outfile`: The path to the output file (must be a gzipped file).

    Example
    -------
    Running the script with:
    
        python script.py --infile data.raw.gz --outfile processed.gz

    Will set:
    - `infile` to `"data.raw.gz"`
    - `outfile` to `"processed.gz"`
    """

    parser = argparse.ArgumentParser(description= "extracting information from raw files for gwas-norm")
    parser.add_argument('--infile', type = str, help = 'the file name of the input file including the ".gz" extension')
    parser.add_argument('--bimfile', type = str, help = 'the file name of the bim file to get the rsids')
    parser.add_argument('--outfile', type = str, help = 'the file name of the output file including the ".gz" extension')

    args = parser.parse_args()
    
    return args.infile, args.bimfile, args.outfile

def format_file(infile, bimfile, outfile, use_cols=['chr_name', 'start_pos', 'effect_allele', 'other_allele', 'uni_id', 'pvalue']):
    """
    Process and format a GWAS summary statistics file for use with magma.
    """
    # Read the bimfile in chunks for memory efficiency
    bimfile_chunksize = 100000  # Adjust chunk size for your memory
    bimfile_iter = pd.read_csv(bimfile, sep='\t', dtype='str', 
                               names=['chr_name', 'rsid', 'start_pos', 'other_allele', 'effect_allele'],
                               usecols=[0, 1, 3, 4, 5], chunksize=bimfile_chunksize)

    # Open a temporary uncompressed file for intermediate results
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_path = temp_file.name

    with gzip.open(infile, 'rt') as infile, open(temp_path, 'w') as outfile_temp:
        reader = pd.read_csv(infile, usecols=use_cols, sep='\t', chunksize=100000)  # Read chunks of the GWAS file
        
        for chunk in reader:
            # Merge with `bimfile` (only necessary columns)
            for bim_chunk in bimfile_iter:
                # Create the unique identifiers for merging (uni_id1 and uni_id2)
                bim_chunk['uni_id1'] = bim_chunk['chr_name'] + '_' + bim_chunk['start_pos'] + '_' + bim_chunk['other_allele'] + '_' + bim_chunk['effect_allele']
                bim_chunk['uni_id2'] = bim_chunk['chr_name'] + '_' + bim_chunk['start_pos'] + '_' + bim_chunk['effect_allele'] + '_' + bim_chunk['other_allele']
                bim_chunk = bim_chunk[['rsid', 'uni_id1', 'uni_id2']]

                # Merge with both `uni_id1` and `uni_id2`
                merge1 = chunk.merge(bim_chunk[['rsid', 'uni_id1']], left_on='uni_id', right_on='uni_id1', how='left')
                merge2 = chunk.merge(bim_chunk[['rsid', 'uni_id2']], left_on='uni_id', right_on='uni_id2', how='left')

                # Concatenate and drop duplicates
                chunk = pd.concat([merge1, merge2]).drop_duplicates().reset_index(drop=True)

                # Apply any necessary transformations
                chunk['pvalue'] = np.exp(-(chunk['pvalue']))

                # Write the processed chunk to the temporary file
                chunk[['rsid', 'pvalue','chr_name','start_pos', 'effect_allele', 'other_allele']].to_csv(outfile_temp, sep='\t', index=False, header=outfile_temp.tell() == 0, mode='a')

                del chunk  # Free memory after each chunk is processed
                break  # Only process one chunk of bimfile at a time to save memory

    # Now compress the temporary uncompressed file into the final output file
    with open(temp_path, 'rb') as f_in, gzip.open(outfile, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Remove the temporary uncompressed file
    os.remove(temp_path)


def main():
    infile, bimfile, outfile = get_args()
    format_file(infile, bimfile, outfile)

if __name__ == '__main__':
    main()