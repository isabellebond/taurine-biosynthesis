import pandas as pd
import os
import gzip
import shutil
import tempfile
import argparse

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
    parser.add_argument('--outfile', type = str, help = 'the file name of the output file including the ".gz" extension')

    args = parser.parse_args()
    
    return args.infile, args.outfile

def format_file(infile, outfile):
    """
    Process and format a GWAS summary statistics file for normalization with gwas-norm.

    This function reads a gzipped input file in chunks, processes the data by extracting
    genomic information from the 'MarkerName' column, optimizes data types for memory efficiency,
    and writes the processed data to a temporary uncompressed file before compressing it again.

    Parameters
    ----------
    infile : str
        Path to the input GWAS summary statistics file (gzipped, ".gz").
    outfile : str
        Path to the output formatted GWAS file (gzipped, ".gz").

    Processing Steps
    ----------------
    - Reads the input file in chunks (100,000 rows at a time).
    - Splits the 'MarkerName' column into 'chromosome', 'start_pos', 'other_allele', and 'effect_allele'.
    - Ensures 'end_pos' is identical to 'start_pos'.
    - Removes the 'chr' prefix from chromosome values.
    - Converts columns to optimized datatypes to reduce memory usage.
    - Writes processed data to a temporary uncompressed file.
    - Compresses the temporary file into the final gzipped output.

    Notes
    -----
    - Uses `tempfile.NamedTemporaryFile` to create a temporary file for intermediate storage.
    - Uses `shutil.copyfileobj` for efficient file compression.
    - Deletes the temporary file after writing the final gzipped output.

    Example
    -------
    >>> format_file("gwas_input.gz", "gwas_output.gz")

    This will process "gwas_input.gz" and save the formatted output as "gwas_output.gz".
    """
     
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_path = temp_file.name

    # Process in chunks and write to an uncompressed temp file
    with gzip.open(infile, 'rt') as infile, open(temp_path, 'w') as outfile_temp:
        reader = pd.read_csv(infile, sep='\t', chunksize=100000)

        for chunk in reader:
            split_cols = chunk['MarkerName'].str.split(':', expand=True)
            chunk[['chromosome', 'start_pos', 'other_allele', 'effect_allele']] = split_cols

            chunk['end_pos'] = chunk['start_pos']  # Duplicate start_pos for end_pos

            # Remove 'chr' prefix and optimize datatypes
            chunk['chromosome'] = chunk['chromosome'].str.replace('chr', '', regex=False)
            chunk = chunk.astype({
                'chromosome': 'category',
                'start_pos': 'int32',
                'end_pos': 'int32',
                'other_allele': 'category',
                'effect_allele': 'category'
            })

            # Write processed chunk to temporary file
            chunk.to_csv(outfile_temp, sep='\t', index=False, header=outfile_temp.tell()==0, mode='a')
            del chunk  # Free memory

    with open(temp_path, 'rb') as f_in, gzip.open(outfile, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Remove temporary uncompressed file
    os.remove(temp_path)


def main():
    infile, outfile = get_args()
    format_file(infile, outfile)

if __name__ == '__main__':
    main()
       
        