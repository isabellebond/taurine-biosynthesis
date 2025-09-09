import pandas as pd
import logging
import scanpy as sc
import gseapy as gp

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def load_counts(file_path, id_col = 'gene_id'):
    """Load gene expression counts from a CSV file."""
    df = pd.read_csv(file_path, index_col=0, sep = '\t')

    genemap = dict(zip(df.index, df[id_col]))

    df = df.drop(columns=[id_col])
    df = df.astype(int)

    logging.info(f"Loaded counts data with {len(df)} genes and {len(df.columns)} samples.")
    
    logging.info(f"Filtered counts data to {len(df)} genes with non-zero counts.")

    df = df.T

    metadata = pd.DataFrame(index=df.index)
    metadata['Differentiation'] = metadata.index.str.split('_').str[0]
    metadata['Genotype'] = metadata.index.str.split('_').str[1]
    metadata['Drug'] = metadata.index.str.split('_').str[2]
    metadata['Combination'] = metadata['Genotype'] + '_' + metadata['Drug']

    print(metadata)

    return df, metadata, genemap

def run_deseq2(counts, metadata, genemap, design_formula = ['Genotype', 'WT', 'TTN'], lowly_expressed = 10, outfile = None):
    #filter out genes with zero counts across all samples
    counts = counts.loc[counts.sum(axis=1) > 0]

    dds = DeseqDataSet(counts = counts,
                   metadata = metadata,
                   design = f'~{design_formula[0]}')
    
    dds.deseq2()

    stat_res = DeseqStats(dds, n_cpus=1, contrast = design_formula)
    stat_res.summary()
    res = stat_res.results_df
    logging.info(f"Differential expression analysis complete for {len(res)} genes.")

    res = res[res.baseMean >= lowly_expressed]

    res['gene_name'] = res.index.map(genemap)
    if outfile:
        res.to_csv(outfile, sep='\t', index=False)
        logging.info(f"Results saved to {outfile}.")

    logging.info(f"{len(res)} genes with baseMean >= 10.")

    sig = res[(res['padj'] < 0.05) & (abs(res['log2FoldChange']) > 0.5)]
    logging.info(f"{len(sig)} significantly differentially expressed genes (padj < 0.05 and |log2FC| > 0.5).")
    top_genes = "\n".join(sig.sort_values("padj")["gene_name"].head(10).tolist())
    logging.info(f"The top 10 significant genes are:\n{top_genes}")

    return res, sig

def run_gsea(res, outfile = None):
    res = res.reset_index()
    ranking = res[['gene_id', 'stat']].dropna().sort_values(by='stat', ascending=False)
    #print duplicates
    duplicates = ranking[ranking.duplicated(subset='gene_id', keep=False)]
    if len(duplicates) > 0:
        logging.warning(f"Found {len(duplicates)} duplicate gene names in ranking. Keeping the first occurrence.")
        print(duplicates)
    ranking = ranking.drop_duplicates(subset='gene_id', keep = 'first')
    logging.info(f"Removed duplicates, {len(ranking)} unique genes remain for GSEA.")

    pre_res = gp.prerank(rnk=ranking, 
                         gene_sets=['GO_Molecular_Function_2025',
                                    'GO_Cellular_Component_2025',
                                    'GO_Biological_Process_2025'], outdir=None, permutation_num=1000, seed=6, verbose=True)
    
    
    terms = pre_res.res2d
    print(terms)
    terms['Dataset'] = terms['Term'].str.split('__').str[0]
    terms['GO_ID'] = terms['Term'].str.split('(').str[1].str.split(')').str[0]
    terms['Term'] = terms['Term'].str.split('__').str[1].str.split('(').str[0].str.strip()

    terms.sort_values(by='FWER p-val', ascending=True, inplace=True)

    if outfile:
        terms.to_csv(outfile, sep='\t', index=False)
        logging.info(f"GSEA results saved to {outfile}.")

    logging.info(f"GSEA complete with {len(sig)} enriched terms.")
    top_terms = "\n".join(sig["Term"].head(10).tolist())
    logging.info(f"The top 10 enriched terms are:\n{top_terms}\n")

    return sig


def main():
    logging.basicConfig(
        filename="logging/deg_rnaseq.log",   # log file name
        filemode="w",                   # overwrite each run ("a" to append)
        format="%(message)s",           # just the message, no timestamps/levels
        level=logging.INFO
    )

    gene_counts, gene_metadata, genemap = load_counts("rnaseq/taurine_alignment_salmon/salmon.merged.gene_counts.tsv", id_col = 'gene_name')
    #transcript_counts, transcript_metadata, transmap = load_counts("rnaseq/taurine_alignment_salmon/salmon.merged.transcript_counts.tsv", id_col = 'gene_id')

    combinations = [['WT_Untreated', 'TTN_Untreated'],
                    ['WT_Untreated', 'TTN_Taurine'],
                    ['WT_Untreated', 'WT_Taurine'], 
                    ['TTN_Untreated', 'TTN_Taurine'], 
                    ['WT_Taurine', 'TTN_Untreated'],
                    ['WT_Taurine', 'TTN_Taurine']]
    
    for combination in combinations:
        logging.info(f"-----------------------------------------------------------")
        logging.info(f"Analyzing combination: {combination[0]} vs {combination[1]}")
        logging.info(f"-----------------------------------------------------------")

        comb_counts = gene_counts[gene_metadata['Combination'].isin(combination)]
        comb_metadata = gene_metadata[gene_metadata['Combination'].isin(combination)]

        res, sig = run_deseq2(comb_counts, comb_metadata, genemap, design_formula = ['Combination', combination[0], combination[1]] , outfile = f'results/deseq2_{combination[0]}_vs_{combination[1]}.tsv')
        run_gsea(res, outfile = f'results/gsea_{combination[0]}_vs_{combination[1]}.tsv')

        #comb_counts = transcript_counts[transcript_metadata['Combination'].isin(combination)]
        #comb_metadata = transcript_metadata[transcript_metadata['Combination'].isin(combination)]
        #res, sig = run_deseq2(comb_counts, comb_metadata, transmap, design_formula = ['Combination', combination[0], combination[1]] , outfile = f'results/deseq2_transcript_{combination[0]}_vs_{combination[1]}.tsv')
        logging.info(f'\n')

    
if __name__ == "__main__":
    main()






        
    