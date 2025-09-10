import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import scanpy as sc
import gseapy as gp
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def load_counts(file_path, id_col = 'gene_id', drop_samples = []):
    """Load gene expression counts from a CSV file."""
    df = pd.read_csv(file_path, index_col=0, sep = '\t')

    genemap = dict(zip(df.index, df[id_col]))

    df = df.drop(columns=[id_col])
    df = df.drop(columns=drop_samples)
    
    df = df.astype(int)

    logging.info(f"Loaded counts data with {len(df)} genes and {len(df.columns)} samples.")

    df = df.T
    print(df)
    

    genes_to_keep = df.columns[df.sum(axis=0) >= 10]
    df = df[genes_to_keep]

    logging.info(f"Filtered counts data to {len(df.columns)} genes with greater than 10 total counts across all samples.")

    metadata = pd.DataFrame(index=df.index)
    metadata['Differentiation'] = metadata.index.str.split('_').str[0]
    metadata['Genotype'] = metadata.index.str.split('_').str[1]
    metadata['Drug'] = metadata.index.str.split('_').str[2]
    metadata['Combination'] = metadata['Genotype'] + '_' + metadata['Drug']


    return df, metadata, genemap

def create_deseq2_object(counts, metadata, design = 'Genotype'):
    dds = DeseqDataSet(counts = counts,
                   metadata = metadata,
                   design = f'~{design}',
                   refit_cooks = True)
    
    dds.deseq2()

    return dds

def de_statistics(dds, contrast, alpha = 0.05):
    ds = DeseqStats(
        dds,
        contrast = contrast,
        alpha = alpha,
        cooks_filter = True,
        independent_filter = True
    )

    ds.summary()

    summary = ds.results_df

    logging.info(f'Differential expression analysis run on {len(summary)} genes.')
    summary = summary.loc[summary['baseMean'] != 0]
    logging.info(f'{len(summary)} genes remain after filtering out zero baseMean.')
    logging.info(f'{len(summary.dropna())} genes have non-NA adjusted p-values (NA p values due to independent filter).')
    logging.info(f'{len(summary.loc[summary["padj"]<0.1])} genes have significant changes in differential expression.')
    logging.info(f'Of these {len(summary.loc[((summary["padj"]<0.05) & (summary["log2FoldChange"]>0))])} genes are upregulated and {len(summary.loc[((summary["padj"]<0.05) & (summary["log2FoldChange"]<0))])} are downregulated.')
    
    summary = summary.sort_values(by='padj', ascending=True)

    return summary

def plot_pca(dds, color_by = 'Combination', outfile_prefix = 'results/pca_plot'):
    sc.tl.pca(dds)
    sc.pl.pca(dds, color = color_by, size = 200, title = 'PCA plot', return_fig = True, show = True)
    plt.savefig(f'{outfile_prefix}.png', bbox_inches='tight')
    return

def plot_volcano(df, labels, label_col=None, title=None, hline=1.3, vline=0, savepath=None):
    df = df.dropna(subset=['padj', 'log2FoldChange'])
    df['-log10(padj)'] = -1 * np.log10(df['padj'])
    fig, ax = plt.subplots(figsize=(6, 10))
    
    # Scatter all points
    ax.scatter(df['log2FoldChange'], df['-log10(padj)'], alpha=0.5)
    
    texts = []
    for label in labels:
        row = df[df[label_col] == label]
        if not row.empty:
            x = row['log2FoldChange'].values[0]
            y = row['-log10(padj)'].values[0]
            ax.scatter(x, y, color='red')
            # Add text with initial position
            texts.append(ax.text(x, y, label, fontsize=9, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=0.01)))

    # Adjust text positions to avoid overlap
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5))

    # Reference lines
    if hline:
        ax.axhline(y=hline, color='grey', linestyle='--')
    if vline:
        ax.axvline(x=vline, color='grey', linestyle='--')

    # Y-axis as p-values
    yticks = ax.get_yticks()
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{10**(-y):.1e}" if y >= 0 else "" for y in yticks])

    # X-axis as fold change
    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{2**x:.2f}" for x in xticks])

    ax.set_xlabel("Fold Change")
    ax.set_ylabel("Adjusted p-value")
    if title:
        ax.set_title(title)

    if savepath:
        plt.savefig(savepath, bbox_inches='tight')
    plt.show()
  
def plot_heatmap(dds, genes, genemap, samples = None, prefix = None, row_order = None):
    dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

    reverse_genemap = {v: k for k, v in genemap.items()}
    ensembl_genes = [reverse_genemap[g] for g in genes if g in reverse_genemap]
    
    dds_heatmap = dds[:, ensembl_genes]
    grapher = pd.DataFrame(dds_heatmap.layers['log1p'].T,
                       index=dds_heatmap.var_names, columns=dds_heatmap.obs_names)
    if row_order:
        grapher = grapher[row_order]

    grapher.index = grapher.index.map(genemap)
    if samples:
        grapher = grapher[samples]

    sns.clustermap(grapher, z_score = 0, cmap = 'RdYlBu_r')
    plt.savefig(f'results/{prefix}_heatmap.png', bbox_inches='tight')

def run_gsea(res, outfile_prefix = None):
    res = res.reset_index()

    ranking = res[['gene_name', 'stat']].dropna().sort_values(by='stat', ascending=False)
    #print duplicates
    duplicates = ranking[ranking.duplicated(subset='gene_name', keep=False)]
    if len(duplicates) > 0:
        logging.warning(f"Found {len(duplicates)} duplicate gene names in ranking. Keeping the first occurrence.")
        print(duplicates)
    ranking = ranking.drop_duplicates(subset='gene_name', keep = 'first')
    logging.info(f"Removed duplicates, {len(ranking)} unique genes remain for GSEA.")

    for gene_set in ['GO_Molecular_Function_2025', 'GO_Cellular_Component_2025', 'GO_Biological_Process_2025', 'Reactome_Pathways_2024']:
        pre_res = gp.prerank(rnk=ranking, 
                             gene_sets=gene_set, 
                             permutation_num=1000, seed=6, verbose=True)

        terms = pre_res.res2d
        terms['Dataset'] = terms['Term'].str.split('__').str[0]
        terms['GO_ID'] = terms['Term'].str.split('(').str[1].str.split(')').str[0]
        terms['Term'] = terms['Term'].str.split('__').str[1].str.split('(').str[0].str.strip()

        terms.sort_values(by='FWER p-val', ascending=True, inplace=True)

        if outfile_prefix:
            terms.to_csv(f'{outfile_prefix}_gene_set.txt', sep='\t', index=False)

    return terms


def main():
    logging.basicConfig(
        filename="logging/deg_rnaseq.log",   # log file name
        filemode="w",                   # overwrite each run ("a" to append)
        format="%(message)s",           # just the message, no timestamps/levels
        level=logging.INFO
    )

    gene_counts, gene_metadata, genemap = load_counts("rnaseq/counts/salmon.merged.gene_counts.tsv", id_col = 'gene_name')
    dds = create_deseq2_object(gene_counts, gene_metadata, design = 'Combination')
    plot_pca(dds, color_by = 'Combination', outfile_prefix='results/pca_plot_all_samples.png')

    gene_counts, gene_metadata, genemap = load_counts("rnaseq/counts/salmon.merged.gene_counts.tsv", id_col = 'gene_name', drop_samples=['431_WT_Taurine'])
    dds = create_deseq2_object(gene_counts, gene_metadata, design = 'Combination')
    plot_pca(dds, color_by = 'Combination', outfile_prefix='results/pca_plot_filtered_samples.png')
    #print(gene_counts)

    combinations = [['WT_Untreated', 'TTN_Untreated'],
                    ['WT_Untreated', 'TTN_Taurine'],
                    ['WT_Untreated', 'WT_Taurine'], 
                    ['TTN_Untreated', 'TTN_Taurine'], 
                    ['WT_Taurine', 'TTN_Taurine']]

    for combination in combinations:
        summary = de_statistics(dds, ['Combination', combination[0], combination[1]])
        summary['gene_name'] = summary.index.map(genemap)

        summary.reset_index(inplace=True)
        summary['gene_name'] = summary['gene_name'].fillna(summary['gene_id'])
        summary.to_csv(f'results/deseq2_{combination[0]}_vs_{combination[1]}.tsv', sep = '\t', index = False)

        #summary = pd.read_csv(f'results/deseq2_{combination[0]}_vs_{combination[1]}.tsv', sep = '\t')
        labels = summary.loc[summary['padj']<0.05, 'gene_name'].tolist()
        
        mitochondrial_genes = summary.loc[summary['gene_name'].str.startswith('MT-'), 'gene_name'].tolist()
        genesets = [labels, 
                    mitochondrial_genes,
                    ['TTN', 'NPPA','NPPB','NDUFC2', 'NDUFS5', 'COX5B', 'MYL2', 'VEGF']
                    ]
        geneset_names = ['sig_deg', 'mitochondrial', 'ttn_important']
        for geneset, name in zip(genesets, geneset_names):
            plot_heatmap(dds, genes = geneset, genemap = genemap, prefix = f'{name}_{combination[0]}_vs_{combination[1]}')
        
        #run_gsea(summary, outfile_prefix = 'results/gsea_WT_vs_TTN_')

    #

    #plot_volcano(summary, labels = labels, label_col = 'gene_name', title = 'WT_Untreated vs TTN_Untreated', 
    #              hline = 1.3, vline = 0, savepath = 'results/volcano_WT_Untreated_vs_TTN_Untreated.png')
    
    #plot_heatmap(summary, rows = labels)

    
    
    #transcript_counts, transcript_metadata, transmap = load_counts("rnaseq/taurine_alignment_salmon/salmon.merged.transcript_counts.tsv", id_col = 'gene_id')
    """
    combinations = [['WT_Untreated', 'TTN_Untreated'],
                    ['WT_Untreated', 'TTN_Taurine'],
                    ['WT_Untreated', 'WT_Taurine'], 
                    ['TTN_Untreated', 'TTN_Taurine'], 
                    ['WT_Taurine', 'TTN_Taurine']]
    
    for combination in combinations:
        logging.info(f"-----------------------------------------------------------")
        logging.info(f"Analyzing combination: {combination[0]} vs {combination[1]}")
        logging.info(f"-----------------------------------------------------------")

        comb_counts = gene_counts[gene_metadata['Combination'].isin(combination)]
        comb_metadata = gene_metadata[gene_metadata['Combination'].isin(combination)]

        res, sig = run_deseq2(comb_counts, comb_metadata, genemap, design = 'Combination' , outfile = f'results/deseq2_{combination[0]}_vs_{combination[1]}.tsv')
        run_gsea(res, outfile = f'results/gsea_{combination[0]}_vs_{combination[1]}.tsv')

        
        if 'WT_Taurine' in combination:
            logging.info(f'Analysis rerun with WT_439_Taurin removed due to outlier status in PCA plot.')
            comb_counts = comb_counts.drop(columns=['WT_439_Taurine'])
            comb_metadata = comb_metadata.drop(index=['WT_439_Taurine'])

        #comb_counts = transcript_counts[transcript_metadata['Combination'].isin(combination)]
        #comb_metadata = transcript_metadata[transcript_metadata['Combination'].isin(combination)]
        #res, sig = run_deseq2(comb_counts, comb_metadata, transmap, design_formula = ['Combination', combination[0], combination[1]] , outfile = f'results/deseq2_transcript_{combination[0]}_vs_{combination[1]}.tsv')
        logging.info(f'\n')
    """
    
if __name__ == "__main__":
    main()






        
    