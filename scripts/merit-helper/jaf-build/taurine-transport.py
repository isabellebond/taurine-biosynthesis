#!/usr/bin/env python
'''
Preparing a small cluster job to evaluate the effects of taurine and its products on DCM and associated traits.

Question 1: Which of taurine or its products is associated with risk of DCM?
a. Perform genome wide MR of taurine + products with DCM


Question 2: For associated products, how is the product acting?
a. Perform cis-MR with eQTLs as exposures and taurine + products as outcome
b. Perform cis-MR with cis-metabolites as exposures and disease as outcome
c. Perform cis-MR with eQTLs as exposures and disease as outcome

'''
###############################################################################
# preparing a job array file
###############################################################################

# imports
import os
import pandas as pd
from pathlib import Path

# constants
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
print(ROOT_DIR)
HOME = str(Path.home())
# NOTE: the output directory
OUTPUT_PATH = os.path.join(HOME, 'results/merit-helper/taurine-transport')
Path(OUTPUT_PATH).mkdir(parents=True, exist_ok=True)

# ordered columns
selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

metabolites = ['taurine', 'taurochlolate', 'acetylcarnitine', 'n-acetyltaurine']
genes = ['PTER','BAAT', 'SLC6A6', 'SLC6A11', 'SLC5A13']

###############################################################################
#Get genome wide inputs
###############################################################################
selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

exposures = pd.read_csv('data/exposure-files/metabolite_exposures.txt', sep = '\t')
metabolite_exposures = exposures.loc[exposures['phenotype'].isin(metabolites)]

jaf = metabolite_exposures.copy()
jaf['infile'] = os.path.join(ROOT_DIR, 'data', 'outcome-files', 'outcomes.txt')
jaf['exposure'] = jaf['phenotype']
jaf['ensemblid'] = "genomewide"
jaf['exposure_path'] = os.path.join(ROOT_DIR, 'data', 'exposure-files', 'metabolite_exposures.txt')
jaf['up'] = 25000
jaf['down'] = 25000
jaf['pvalue'] = 8
jaf['logpval'] = True
jaf['maf'] = 0.05
jaf['ld_cut'] = '0.05'
jaf['ld_sample'] = 10000
jaf['ld_seed'] = 12052018
jaf['ref_pop'] = 'EUR_UKB_GRCh38_without_homozygosity'
jaf['sample_list'] = 'EUR_UKB_WO_RELATED'
jaf['proximity_dist'] = 'None'
jaf['drop_variants'] = 'None'
jaf["steiger_pvalue"] = 0

# #### Which models to use
# will perform an IVW, IVW pruned IVW, Egger, Egger pruned Egger
jaf['models'] = "IVW;IVW|IVW;Egger;Egger|Egger"
jaf['models-kwargs'] = 'None'
jaf['outfile'] = [OUTPUT_PATH + '/' +  str(i).strip() + '.tar.gz' for\
                            i in jaf['phenotype'].to_list()]

jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))

jaf = jaf[selcol]
old_jaf = jaf.copy()


###############################################################################
#Run cis-mrs with the same metabolite data
###############################################################################

#Get ensembl info for genes of interest
ensembl = pd.read_csv('data/merit-helper/mapping-files/ensembl.genes.txt', sep = '\t')
ensembl = ensembl.loc[ensembl['gene_name'].isin(genes)]

jafs = []

for i,row in ensembl.iterrows():
    gene_jaf = old_jaf.copy()
    gene_jaf['pvalue'] = 6
    gene_jaf['ensemblid'] = row['gene_id']
    gene_jaf['phenotype'] = gene_jaf['exposure'] + '_' + row['gene_name'] 
    gene_jaf['outfile'] = [OUTPUT_PATH + '/' +  str(i).strip() + '.tar.gz' for\
                            i in jaf['phenotype'].to_list()]
    jaf = pd.concat([jaf, gene_jaf])

jaf.to_csv('data/merit-helper/jaf-files/taurine-transport.jaf', sep = '\t', index = False)