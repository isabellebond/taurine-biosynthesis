import pandas as pd
import warnings
import scanpy as sc

H5AD_FILE = "data/processed/scRNAseq/all_samples_filtered.h5ad"
# Donor IDs with TTN mutations
TTN_DONOR_IDS = [1290, 1304, 1371, 1472] 
# Donor IDs with no known CM mutations
NF_DONOR_IDS = [1515, 1516, 1539, 1540, 1547,
                1558, 1561, 1610, 1678, 1702,
                1718] #1603,1606,1622 removed following QC (ST2)
CELL_TYPE = ['Cardiomyocyte_I', 'Cardiomyocyte_II', 'Cardiomyocyte_III']
CELL_GROUP = 'Cardiomyocyte'
UMI_THRESHOLD = 150
GENE_THRESHOLD = 150
MITO_THRESHOLD = 5
CELL_THRESHOLD = 25

def load_and_filter_data(file_path):
    adata = sc.read_h5ad(file_path)
    
    # Filter for cardiomyocytes
    adata = adata[adata.obs['cell_type_leiden0.6'].isin(CELL_TYPE)]
    adata = adata[adata.obs['donor_id'].isin(TTN_DONOR_IDS + NF_DONOR_IDS)]
    adata = adata[adata.obs['cellranger_percent_mito'] < 5]
    adata = adata[adata.obs['cellbender_ncount'] > UMI_THRESHOLD]
    adata = adata[adata.obs['cellbender_ngenes'] > GENE_THRESHOLD]

    #group cardiomyocytes into one cell type
    adata.obs['cell_type'] = CELL_GROUP

    return adata

def pseudobulk_aggregation(adata, layer = None):

    #pseudobulk filtering
    cell_counts = adata.obs['donor_id'].value_counts()
    valid_donors = cell_counts[cell_counts >= CELL_THRESHOLD].index

    #psuedobulk aggregation for cellbender counts
    pbs = []
    for donor in adata.obs.donor_id.unique():
        subset = adata[(adata.obs['donor_id'] == donor)]
        if subset.shape[0] < CELL_THRESHOLD:
            warnings.warn(f"Skipping donor {donor} with only {subset.shape[0]} cells")
            continue

        if layer:
            rep_adata = sc.AnnData(X=subset.layers[layer].sum(axis=0),
                                var=subset.var)
        else:
            rep_adata = sc.AnnData(X=subset.X.sum(axis=0),
                                var=subset.var)
        pbs.append(rep_adata)
    
    pbs = sc.concat(pbs)
    counts = pd.DataFrame(pbs.X, columns=adata.var_names)
    meta = rep_adata.obs
    
    return counts, meta

def main():
    adata = load_and_filter_data(H5AD_FILE)

    counts, meta = pseudobulk_aggregation(adata)
    counts.to_csv("results/pseudobulk_counts.cellbender.csv")
    meta.to_csv("results/pseudobulk_metadata.csv")

    counts, meta = pseudobulk_aggregation(adata, layer='cellranger_raw')
    counts.to_csv("results/pseudobulk_counts.cellranger.csv")

if __name__ == "__main__":
    main()

