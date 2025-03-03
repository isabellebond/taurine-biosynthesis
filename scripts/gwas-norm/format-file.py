import pandas as pd
import os
import gzip
import shutil
import tempfile
import argparse

def get_args():
    parser = argparse.ArgumentParser(description= "extracting information from raw files for gwas-norm")
    parser.add_argument('--infile', type = str, help = 'the file name of the input file including the ".gz" extension')
    parser.add_argument('--outfile', type = str, help = 'the file name of the output file including the ".gz" extension')

    args = parser.parse_args()
    
    return args.infile, args.outfile

def format_file(infile, outfile):
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
       
        