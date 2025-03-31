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
ROOT_DIR = os.getcwd()
print(ROOT_DIR)
HOME = str(Path.home())
# NOTE: the output directory
OUTPUT_PATH = os.path.join(HOME, 'Scratch/results_dir/taurine-biosynthesis')
Path(OUTPUT_PATH).mkdir(parents=True, exist_ok=True)

# ordered columns
selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

metabolites = ['Cystine', 'Cysteine', 'Cysteine Sulfinic Acid', 'Taurine', 'Hypotaurine']
genes = ['GCLC','CDO1','GAD1', 'GADL1', 'CSAD', 'CTH','ADO', 'VNN1', 'VNN2', 'VNN3',
         'FMO1','FMO2','FMO3','FMO4', 'FMO5', 'HTDH']

###############################################################################
#Get genome wide inputs
###############################################################################
selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

exposures = pd.read_csv('/home/rmgpibo/taurine-biosynthesis/data/gwas-norm/metabolite-map.txt', sep = '\t')
exposures.rename(columns = {'Path':'path',
                            'Metabolite':'phenotype',
                            'Samplesize':'samplesize'}, inplace = True)
exoposures = exposures.loc[exposures['phenotype'].isin(metabolites)]
exposures['phenotype'] = exposures['phenotype'].str.replace(' ','_')

exposure_cols = ['path','effect_type','effect_allele','other_allele','effect_size','standard_error',
 	'pvalue','chr_name','start_pos','end_pos','var_id','chrpos','chrpos_spec','start_anchor',
    'pvalue_logged','sep','compression','group','drop','notes','phenotype','unit','no_cases',
    'samplesize','url']

exposures['effect_type'] = 'beta'
exposures['effect_allele'] = None
exposures['other_allele'] = None
exposures['effect_size'] = None
exposures['standard_error'] = None
exposures['pvalue'] = None
exposures['chr_name'] = None
exposures['start_pos'] = None
exposures['end_pos'] = None
exposures['chrpos'] = None
exposures['chrpos_spec'] = None
exposures['start_anchor'] = None
exposures['end_anchor'] = None
exposures['pvalue_logged'] = True
exposures['sep'] = '\t'
exposures['compression'] = 'gzip'
exposures['group'] = None
exposures['drop'] = None
exposures['notes'] = None
exposures['var_id'] = None
exposures['unit'] = None
exposures['no_cases'] = 0
exposures['url'] = None

expsures = exposures[exposure_cols]
print(exposures)

jaf = exposures.copy()
jaf['infile'] = os.path.join(ROOT_DIR, 'data', 'merit-helper','outcome-files', 'outcomes.txt')
jaf['exposure'] = jaf['phenotype']
jaf['ensemblid'] = "genomewide"
jaf['exposure_path'] = os.path.join(ROOT_DIR, 'data', 'merit-helper','exposure-files', 'taurine_biosynthesis.txt')
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
jaf.to_csv('data/merit-helper/jaf-files/taurine-gw.jaf', sep = '\t', index = False)
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
                            i in gene_jaf['phenotype'].to_list()]
    jafs.append(gene_jaf)

jaf_cis = pd.concat(jafs)

###############################################################################
#Run cis-mrs with the eqtl data for given genes
###############################################################################
#Get GTEx data for genes of interest
gtex =  pd.read_csv('data/merit-helper/mapping-files/eqtlgen_eqtl_mapper.tsv', sep = '\t')
eqtlgen = pd.read_csv('data/merit-helper/mapping-files/gtex_eqtl_mapper.tsv', sep = '\t')

gtex = gtex.loc[gtex['ensembl_id'].isin(ensembl['gene_id'].unique())]
eqtlgen = eqtlgen.loc[eqtlgen['ensembl_id'].isin(ensembl['gene_id'].unique())]

eqtls = pd.concat([gtex,eqtlgen])
eqtls.set_index('Unnamed: 0', drop = False, inplace = True)
print(eqtls.columns)
for i,row in eqtls.iterrows():
    print(row)
    gene_jaf = old_jaf.head(1).copy()
    gene_jaf['exposure'] = row['Unnamed: 0']
    gene_jaf['exposure_path'] = os.path.join(ROOT_DIR, 'data','merit-helper', 'exposure-files', 'taurine_biosynthesis_eqtl.txt')
    gene_jaf['pvalue'] = 6
    gene_jaf['ensemblid'] = row['ensembl_id']
    gene_jaf['phenotype'] = gene_jaf['exposure']
    gene_jaf['outfile'] = [OUTPUT_PATH + '/' +  str(i).strip() + '.tar.gz' for\
                            i in gene_jaf['phenotype'].to_list()]
    jaf = pd.concat([jaf, gene_jaf])

jaf = pd.concat([jaf_cis,jaf])
jaf = jaf.drop_duplicates(keep = 'first')
jaf = jaf.reset_index(drop = True)
jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))

exposures = exposures[exposure_cols]
eqtls = eqtls[exposure_cols]

###############################################################################
#Run cis-mrs with the pqtl data for given genes
###############################################################################
#Get pqtl data for genes of interest
pqtl =  pd.read_csv('data/merit-helper/mapping-files/dtadb_protein_mapping_b38.2.txt', sep = '\t')

pqtl = pqtl.loc[pqtl['ensembl_gene_id'].isin(ensembl['gene_id'].unique())]

print(pqtl)

for i,row in pqtl.iterrows():
    print(row)
    gene_jaf = old_jaf.head(1).copy()
    gene_jaf['exposure'] = row['gene'] + '_' + row['dataset']
    gene_jaf['exposure_path'] = os.path.join(ROOT_DIR, 'data','merit-helper', 'exposure-files', 'taurine_biosynthesis_pqtl.txt')
    gene_jaf['pvalue'] = 6
    gene_jaf['ensemblid'] = row['ensembl_gene_id']
    gene_jaf['phenotype'] = gene_jaf['exposure']
    gene_jaf['outfile'] = OUTPUT_PATH + '/' +  row['gene'] + '_' + row['dataset']  + '.tar.gz'
    jaf = pd.concat([jaf, gene_jaf])

pqtl['effect_type'] = 'beta'
pqtl['effect_allele'] = None
pqtl['other_allele'] = None
pqtl['effect_size'] = None
pqtl['standard_error'] = None
pqtl['pvalue'] = None
pqtl['chr_name'] = None
pqtl['start_pos'] = None
pqtl['end_pos'] = None
pqtl['chrpos'] = None
pqtl['chrpos_spec'] = None
pqtl['start_anchor'] = None
pqtl['end_anchor'] = None
pqtl['pvalue_logged'] = True
pqtl['sep'] = '\t'
pqtl['compression'] = 'gzip'
pqtl['group'] = None
pqtl['drop'] = None
pqtl['notes'] = None
pqtl['var_id'] = None
pqtl['unit'] = None
pqtl['no_cases'] = 0
pqtl['url'] = None

jaf = jaf.drop_duplicates(keep = 'first')
jaf = jaf.reset_index(drop = True)
jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))


exposures = exposures[exposure_cols]
exposures.set_index('phenotype', drop = False, inplace = True)


pqtl['phenotype'] = pqtl['gene'] + '_' + pqtl['dataset']
pqtl = pqtl[exposure_cols]
pqtl.set_index('phenotype', drop = False, inplace = True)

exposures.to_csv('data/merit-helper/exposure-files/taurine_biosynthesis.txt', sep = '\t')
jaf.to_csv('data/merit-helper/jaf-files/taurine-biosynthesis.jaf', sep = '\t', index = False)
eqtls.to_csv(os.path.join('data', 'merit-helper', 'exposure-files', 'taurine_biosynthesis_eqtl.txt'), sep = '\t')
pqtl.to_csv(os.path.join('data', 'merit-helper', 'exposure-files', 'taurine_biosynthesis_pqtl.txt'), sep = '\t')
