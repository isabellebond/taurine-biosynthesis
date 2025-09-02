import pandas as pd 
from single_cell_rnaseq.example import load_h5ad, make_pseudobulk, get_deg
import diffxpy.api as de
import scanpy as sc

def get_mouse_data(gene_df):
    mouse = pd.read_csv('data/mouse_data/mouse_phenotypes_opentargets.csv')
    
    # Split into cardiac vs other
    mouse_cardiac = mouse.loc[mouse['Parent Name'] == 'cardiovascular system phenotype']
    mouse_other = mouse.loc[mouse['Parent Name'] != 'cardiovascular system phenotype']

    # Merge cardiac
    df = gene_df.merge(mouse_cardiac[['Phenotype', 'Ensembl Gene ID']], 
                       on='Ensembl Gene ID', how='left')
    df = df.groupby(gene_df.columns.tolist(), dropna=False).agg({
        'Phenotype': lambda x: ', '.join(x.dropna().unique())
    }).reset_index()
    df = df.rename(columns={'Phenotype': 'Cardiac Phenotype'})
    
    # Merge other
    df = df.merge(mouse_other[['Phenotype', 'Ensembl Gene ID']], 
                  on='Ensembl Gene ID', how='left')
    df = df.groupby(gene_df.columns.tolist() + ['Cardiac Phenotype'], dropna=False).agg({
        'Phenotype': lambda x: ', '.join(x.dropna().unique())
    }).reset_index()
    df = df.rename(columns={'Phenotype': 'Other Phenotype'})
    
    return df

def get_differential_expression_data(gene_df):
    genes = list(gene_df['gene_name'].unique())
    print(genes)
    load_h5ad('/Volumes/Ultra Touch/35732739_rnaseq/anndata/human_dcm_hcm_scportal_03.17.2022.h5ad', disease_subset = None, gene_subset = genes, save_path = f'data/deg/taurine_proteins.h5ad')
    
    return

def format_deg_data():
    adata = load_h5ad('/Volumes/Ultra Touch/35732739_rnaseq/anndata/human_dcm_hcm_scportal_03.17.2022.h5ad', disease_subset=None, )
    print(adata.obs['disease'].unique())
    subset = adata

    subset.X = subset.X.toarray()  # Convert sparse matrix to dense'=

    pseudobulk, logcpm,  meta = make_pseudobulk(subset)
    print(meta.head())
    print(meta['disease'].unique())

    all_deg = []

    for ct in meta.index.get_level_values("cell_type_leiden0.6").unique():
        # subset counts + metadata for one cell type
        counts_ct = logcpm.xs(ct, level="cell_type_leiden0.6")
        meta_ct = meta.xs(ct, level="cell_type_leiden0.6")

        # run DE
        deg_ct = get_deg(counts_ct, meta_ct, use_logcpm=True, disease_ref="NF")
        deg_ct["cell_type"] = ct  # add cell type annotation

        all_deg.append(deg_ct)

    # concatenate results
    deg_results = pd.concat(all_deg, ignore_index=True)
    deg_results.to_csv('data/deg/taurine_proteins_deg.csv', index=False)

    

    """try:
            res = de.test.wald(data = cell_subset, formula_loc='~ 1 + disease', factor_loc_totest='disease')
            res = res.summary()
            res['cell_type'] = celltype
            print(res)
        except Exception as e:
            print(f'Error processing cell type {celltype}: {e}')

        results.append(res)
    results = pd.concat(results, ignore_index=True)
    results.to_csv('data/deg/taurine_proteins_deg.csv', index=False)"""
def main():
    
    #format_deg_data()
    taurine = pd.read_csv('data/taurine_genes.txt', sep='\t')
    get_differential_expression_data(taurine)

    """

    taurine = taurine.rename(columns={'ensembl_id': 'Ensembl Gene ID'})
    taurine = taurine[['Ensembl Gene ID', 'gene_name', 'uniprot_id', 'GO_ID', 'GO_Label']]
    agg_cols = ['GO_ID', 'GO_Label']
    group_cols = [c for c in taurine.columns if c not in agg_cols]

    taurine = taurine.groupby(group_cols, dropna=False).agg({
        'GO_ID': lambda x: ', '.join(x.dropna().unique()),
        'GO_Label': lambda x: ', '.join(x.dropna().unique())
    }).reset_index()

    print(taurine.head())

    taurine = get_mouse_data(taurine)
    taurine = get_differential_expression_data(taurine)

    taurine.to_csv('data/taurine_genes/taurine_mouse_phenotypes.csv', index=False)
    """
if __name__ == "__main__":
    main()

