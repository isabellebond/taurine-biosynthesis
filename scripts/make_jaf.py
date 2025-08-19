import pandas as pd 
import os
from funcs.ontology import load_ontology, search_term

def load_gaf_file(gaf_file):
    gaf_columns = [
        "DB",                   # 1. Database
        "DB_Object_ID",         # 2. Identifier in the database
        "DB_Object_Symbol",     # 3. Gene symbol
        "Relation",            # 4. Qualifier (e.g., NOT, contributes_to)
        "GO_ID",                # 5. GO term ID
        "DB_Reference",         # 6. Reference (e.g., PMID)
        "Evidence_Code",        # 7. Evidence code (e.g., IDA, IEA)
        "With_or_From",         # 8. With/From field
        "Aspect",               # 9. Aspect (F = function, P = process, C = component)
        "DB_Object_Name",       # 10. Full gene name
        "DB_Object_Synonym",    # 11. Synonyms
        "DB_Object_Type",       # 12. Type (e.g., gene, protein)
        "Taxon",                # 13. Taxon ID(s)
        "Date",                 # 14. Date (YYYYMMDD)
        "Assigned_By",          # 15. Annotation provider
        "Annotation_Extension", # 16. Optional annotation extensions
        "Gene_Product_Form_ID"  # 17. Optional gene product isoform ID
    ]
    gaf_df = pd.read_csv(gaf_file, sep='\t', comment='!', header=None, names=gaf_columns)
    return gaf_df

def get_taurine_terms():
    onto, namespace = load_ontology('http://purl.obolibrary.org/obo/go/extensions/go-plus.owl', 'http://purl.obolibrary.org/obo/')
    taurine_terms = search_term(onto, namespace, 'taurine')

    gaf = load_gaf_file('data/gene_ontology/goa_human.gaf')
    gaf = gaf.merge(taurine_terms, on = 'GO_ID', how='inner')
    gaf = gaf[['DB','DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'GO_Label', 'Relation']]
    gaf = gaf.drop_duplicates(keep = 'first')

    ensembl = pd.read_csv('data/merit-helper/mapping-files/ensembl.genes.txt', sep='\t')
    gaf = gaf.merge(ensembl, left_on='DB_Object_ID', right_on='uniprot_id', how='left')
    gaf = gaf.rename(columns={'DB_Object_ID': 'uniprot_id', 'DB_Object_Symbol': 'gene_name', 
                              'gene_id': 'ensembl_id'})
    gaf = gaf[['uniprot_id', 'gene_name', 'ensembl_id', 'entrez_id', 'chromosome', 'start', 'end', 'Relation', 'GO_ID', 'GO_Label']]

    gaf.to_csv('data/gene_ontology/taurine_genes.txt', sep='\t', index=False)

def create_exposure(gene_df, qtls, metabolites):
    exposures = gene_df.copy()

    metabolites_list = [metabolites]
    for ensembl_id in gene_df['ensembl_id'].unique():
        print(ensembl_id)
        metabolite = metabolites.copy()
        metabolite['ensembl_id'] = ensembl_id
        metabolite['phenotype'] = ensembl_id + '_' + metabolite['phenotype']
        metabolites_list.append(metabolite)
    
    exposures = pd.concat([qtls, pd.concat(metabolites_list)], ignore_index=True)
    
    print(exposures)

    exposures['effect_type'] = 'beta'
    exposures['pvalue_logged'] = True
    exposures['compression'] = 'gzip'
    exposures['unit'] = 'std'
    exposures['sep'] = '\t'

    exposure_cols = ['path','effect_type','effect_allele','other_allele','effect_size','standard_error',
 	'pvalue','chr_name','start_pos','end_pos','var_id','chrpos','chrpos_spec','start_anchor',
    'pvalue_logged','sep','compression','group','drop','notes','phenotype','unit','no_cases',
    'samplesize','url']
    
    for col in exposure_cols:
        if col not in exposures.columns:
            exposures[col] = None

    exposures = exposures.drop_duplicates(subset = 'phenotype', keep = 'first')
    exposures = exposures[exposure_cols]
    exposures.index = exposures['phenotype']

    print(exposures)

    return exposures

def make_jaf_mr(exposures, exposurefile, outcomefile, outpath, pvalue=6, r2=0.01, genomewide=False):
    selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

    jaf = exposures.copy()

    jaf['infile'] = outcomefile
    jaf['exposure'] = jaf['phenotype']
    jaf['exposure_path'] = exposurefile

    
    jaf['ensemblid'] = 'genomewide'
    jaf['up'] = 200000
    jaf['down'] = 200000

    for i, row in jaf.iterrows():
        if row['phenotype'].startswith('ENSG'):
            jaf.loc[i, 'ensemblid'] = row['phenotype'].split('_')[0]

            

    jaf['pvalue'] = pvalue
    jaf['logpval'] = True
    jaf['maf'] = 0.01
    jaf['ld_cut'] = str(r2)
    jaf['ld_sample'] = 10000
    jaf['ld_seed'] = 12052018
    jaf['ref_pop'] = 'EUR_UKB_GRCh38_without_homozygosity'
    jaf['sample_list'] = 'EUR_UKB_WO_RELATED'
    jaf['proximity_dist'] = 'None'
    jaf['drop_variants'] = 'None'

    if genomewide:
        jaf['models'] = "IVW;IVW|IVW;Egger;Egger|Egger"
        jaf["steiger_pvalue"] = 0.05
        jaf = jaf.loc[jaf['ensemblid'] == 'genomewide', :]

    else:
        jaf['models'] = "IVW;IVW|IVW"
        jaf["steiger_pvalue"] = 0
        jaf = jaf.loc[jaf['ensemblid'] != 'genomewide', :]

    jaf['models-kwargs'] = 'None'
    
    jaf['drop_variants'] = None

    jaf['outfile'] = [outpath + '/' +  str(i).strip() + '.tar.gz' for\
                            i in jaf['phenotype'].to_list()]

    jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))
    jaf = jaf[selcol]
    jaf = jaf.astype(str)

    return jaf

def make_jaf_coloc(exposures, exposurefile, outcomefile, outpath, pvalue=6, r2=0.01):
    columns = [
    'rowidx', 'infile', 'outfile', 'exposure', 'exposure_path', 'coloc_method',
    'coloc_mode', 'coloc_maxhits', 'coloc_significance', 'coloc_r2_independent',
    'coloc_priors', 'up', 'down', 'pvalue', 'maf', 'ensemblid', 'logpval',
    'ld_sample', 'ld_seed', 'drop_variants', 'ref_pop', 'sample_list'
    ]

    jaf = exposures.copy()
    jaf = jaf.loc[jaf['phenotype'].str.startswith('ENSG'), :]
    jaf['ensemblid'] = jaf['phenotype'].str.split('_').str[0]

    jaf['infile'] = outcomefile
    jaf['exposure'] = jaf['phenotype']
    jaf['exposure_path'] = exposurefile

    # No need for a seperate exposure.txt file
    jaf['coloc_method'] = 'cond;cond'
    jaf['coloc_mode'] = 'allbutone'
    jaf['coloc_maxhits'] = 5
    jaf['coloc_significance'] = pvalue
    jaf['coloc_r2_independent'] = r2
    jaf['coloc_priors'] = '1e-4;1e-4;5e-6'

    jaf['up'] = 200000
    jaf['down'] = 200000

    jaf['pvalue'] = 0
    jaf['logpval'] = True
    jaf['maf'] = 0.01
    jaf['ld_cut'] = '0.05'
    jaf['ld_sample'] = 10000
    jaf['ld_seed'] = 12052018
    jaf['ref_pop'] = 'EUR_UKB_GRCh38_without_homozygosity'
    jaf['sample_list'] = 'EUR_UKB_WO_RELATED'
    jaf['drop_variants'] = None

    jaf['outfile'] = [outpath + '/' +  str(i).strip() + '.tar.gz' for\
                            i in jaf['phenotype'].to_list()]

    jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))
    jaf = jaf[columns]
    jaf = jaf.astype(str)

    return jaf

def make_jaf_cis(df, exposures, exposurefile, outcomefile, outpath, pvalue=6, r2=0.01):
    selcol =['rowidx', 'infile', 'outfile', 'exposure', 'exposure_path',
         'ensemblid',  'up', 'down', 'pvalue', 'logpval', 'maf',
         'ld_cut', 'ld_sample','ld_seed', 'ref_pop',
         'sample_list', 'proximity_dist', 'drop_variants',
         'models', 'models-kwargs', 'steiger_pvalue',
         ]

    jaf = df.copy()

    jaf['infile'] = outcomefile
    jaf['exposure'] = jaf['phenotype']
    jaf['exposure_path'] = exposurefile
    jaf['ensemblid'] = jaf['ensembl_id']

    jaf['up'] = None
    jaf['down'] = None

    for i, row in jaf.iterrows():
        # Check for missing or invalid chromosome
        chr_name = str(row['chr_name_gene']).replace('.0','')

        if chr_name not in chromosome_sizes:
            print(f"Warning: Chromosome {chr_name} not found in chromosome sizes. Skipping {row['ensemblid']}.")
            jaf.drop(i, inplace=True)
            continue

        # Calculate 'up'
        start_pos = row['start_pos_gene']
        jaf.loc[i, 'up'] = str(int(min(200000, start_pos - 1)))

        # Calculate 'down'
        end_pos = row['end_pos_gene']
        chr_end = chromosome_sizes[chr_name]
        jaf.loc[i, 'down'] = str(int(min(200000, chr_end - end_pos - 1)))

    jaf['pvalue'] = pvalue
    jaf['logpval'] = True
    jaf['maf'] = 0.01
    jaf['ld_cut'] = str(r2)
    jaf['ld_sample'] = 10000
    jaf['ld_seed'] = 12052018
    jaf['ref_pop'] = 'EUR_UKB_GRCh38_without_homozygosity'
    jaf['sample_list'] = 'EUR_UKB_WO_RELATED'
    jaf['proximity_dist'] = 'None'
    jaf['drop_variants'] = 'None'
    jaf['models'] = "IVW;IVW|IVW"
    jaf['models-kwargs'] = 'None'
    jaf["steiger_pvalue"] = 0
    jaf['drop_variants'] = None

    jaf['outfile'] = [outpath + '/' +  str(i).strip() + '.tar.gz' for\
                            i in jaf['ensemblid'].to_list()]

    jaf['rowidx'] = list(range(1, jaf.shape[0] + 1))
    jaf = jaf[selcol]
    jaf = jaf.astype(str)

    return jaf

def main():
    #get_taurine_terms()
    gaf = pd.read_csv('data/gene_ontology/taurine_genes.txt', sep='\t')

    all_qtls = []
    for file in os.listdir('data/merit-helper/mapping'):
        dataset = file.split('_')[0]
        qtl = pd.read_csv(f'data/merit-helper/mapping/{file}', sep='\t')

        qtl['dataset'] = dataset
        
        if dataset == 'gtex':
            qtl['phenotype'] = qtl['ensembl_id'] + '_' + qtl['dataset'] + '_' + qtl['tissue']
        else:
            qtl['phenotype'] = qtl['ensembl_id'] + '_' + qtl['dataset']

        qtl = qtl.loc[qtl['ensembl_id'].isin(gaf['ensembl_id']), ['path', 'phenotype', 'ensembl_id', 'samplesize']]

        all_qtls.append(qtl)

    all_qtls = pd.concat(all_qtls, ignore_index=True)


    metabolites = pd.read_csv('data/merit-helper/exposure-files/taurine_biosynthesis.txt', sep='\t')
    metabolites['ensembl_id'] = 'genomewide'
    metabolites = metabolites[['path', 'phenotype','ensembl_id', 'samplesize']]

    
    exposures = create_exposure(gaf, all_qtls, metabolites)
    exposure_path = '/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/exposure-files/exposures.txt'
    exposures.to_csv(exposure_path, sep='\t', encoding="utf-8")

    outcome_path = '/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/outcome-files/outcomes.txt'

    outcomes = pd.read_csv(outcome_path, sep='\t')
    print(outcomes['sep'])

    outpath_mr = '/myriadfs/home/rmgpibo/Scratch/results/taurine/mr'
    outpath_mr_gw = '/myriadfs/home/rmgpibo/Scratch/results/taurine/genomewide_mr'
    outpath_coloc = '/myriadfs/home/rmgpibo/Scratch/results/taurine/coloc'

    jaf = make_jaf_mr(exposures, exposure_path, outcome_path, outpath_mr, pvalue=6, r2=0.01, genomewide=False)
    jaf.to_csv('/myriadfs/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/jaf/cis_mr.jaf', sep='\t', index = False,  encoding="utf-8")

    jaf = make_jaf_mr(exposures, exposure_path, outcome_path, outpath_mr_gw, pvalue=8, r2=0.01, genomewide=True)
    jaf.to_csv('/myriadfs/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/jaf/genomewide_mr.jaf', sep='\t', index = False, encoding="utf-8")

    jaf = make_jaf_coloc(exposures, exposure_path, outcome_path, outpath_coloc, pvalue=6, r2=0.01)
    jaf.to_csv('/myriadfs/home/rmgpibo/Scratch/taurine-biosynthesis/data/merit-helper/jaf/coloc.jaf', sep='\t', index = False, encoding="utf-8")

if __name__ == "__main__":
    main()
    