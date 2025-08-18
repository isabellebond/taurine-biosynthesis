import pandas as pd

from funcs.datahandling import extract_coloc_results

def main():
    """
    coloc = extract_coloc_results(
        folder='/home/rmgpibo/Scratch/results/taurine/coloc/',
        outfile='/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/results/coloc.tsv'
    )

    coloc = coloc.loc[coloc['pp_h4']> 0.8, :]
    coloc['phenotype'] = coloc['gene']
    coloc['gene'] = coloc['gene'].str.split('_').str[0]    
    coloc.set_index('gene', inplace=True)

    ensembl = pd.read_csv(
        '/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/mapping-files/ensembl.genes.txt',
        sep='\t'
    ).set_index('gene_id')

    coloc = coloc.join(ensembl[['gene_name', 'uniprot_id']], how='left')
    coloc = coloc[['gene_name', 'uniprot_id', 'phenotype', 'trait', 'pp_h4']]

    coloc.to_csv(
        '/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/results/coloc.formatted.tsv',
        sep='\t'
    )
    """

    ensembl = pd.read_csv(
        '/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/mapping-files/ensembl.genes.txt',
        sep='\t'
    ).set_index('gene_id')

    mr = pd.read_csv('/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/results/mr.tsv', sep='\t')
    #mr = mr.loc[mr['outcome_crude'].isin(['lvef', 'dcm'])]

    # Count rows per Outcome
    counts = mr.groupby('outcome_crude').size().rename('n')
    print(counts)

    # Merge counts into dataframe
    mr = mr.merge(counts, on='outcome_crude', how = 'left')

    # Calculate per-outcome Bonferroni threshold
    mr['bonf_thresh'] = 0.05 / mr['n']

    # Boolean column for passing Bonferroni
    mr['Bonferroni'] = mr['P-value (float)'] < mr['bonf_thresh']

    mr = mr.loc[mr['Bonferroni'] == True]

    # Optional: if you still want exposure/gene columns
    mr['exposure'] = mr['filename'].str.split('.').str[0]
    mr['ensembl_id'] = mr['exposure'].str.split('_').str[0]
    mr.set_index('ensembl_id', inplace=True)

    mr['outcome'] = mr['outcome_crude']

    mr = mr.join(ensembl[['gene_name', 'uniprot_id']], how='left').reset_index()
    mr.rename(columns={'index': 'ensembl_id'}, inplace=True)
    print(mr)

    mr = mr[['exposure', 'outcome', 'ensembl_id','gene_name', 'uniprot_id', 'No. variants', 'Point estimate', 'Standard error', 'Upper bound', 'Lower bound', 'P-value (float)', 'Bonferroni', 'bonf_thresh', 'n']]
    mr.sort_values(by = ['P-value (float)'], inplace=True, ascending=True)


    mr.to_csv('data/merit-helper/results/mr.formatted.tsv', sep = '\t')

if __name__ == "__main__":
    main()