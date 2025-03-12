#!/bin/bash

run_command() {
    # Decompose the line array, this is really for code clarity
    # Make sure order matches columns in jaf.txt file (ordered in selcol object in jaf_preparation.py)
    INFILE="${PROCESS_LINE[1]}" #GWAS sumstats in line 1
    N="${PROCESS_LINE[2]}" #Number of samples in line 2
    ANNOTFILE="${PROCESS_LINE[3]}"
    OUTFILE="${PROCESS_LINE[4]}" #Output file for analysis

    set +u
    if [ -z "$PYTHON_INTERPRETER" ]; then
         PYTHON_INTERPRETER="${HOME}/miniconda3/envs/taurine/bin/python"
    fi
    set -u

    MUNGE_FEATURES="/lustre/projects/mol_cardio/resources/software/pops/example/data/features_munged"
    POPS="/lustre/projects/mol_cardio/resources/software/pops/pops.py"
    NCBI38="/lustre/projects/mol_cardio/resources/data/magma/NCBI38.gene.loc"
    G1000_EUR="/lustre/projects/mol_cardio/resources/data/magma/g1000_eur.bim"

    # Create SNP location file using awk
    awk 'NR==FNR {a[$2] = $1} NR>FNR {if ($1":"$2":"$4":"$5 in a) print a[$1":"$2":"$4":"$5], $1, $2}' OFS="\t" "${G1000_EUR}.bim" "${INFILE}" > "${ANNOTFILE}_snploc.txt"

    # Create annotation file using magma
    magma --annotate \
            --snp-loc "${TEMP_APP1OUT}_snploc.txt" \
            --gene-loc "${NCBI38}" \
            --out "${ANNOTFILE}"
    
    #get snp and pvalue pairs for analysis
    awk 'NR==FNR {a[$2] = $1} NR==1 {print "SNP", "P"} NR>FNR {if ($1":"$2":"$4":"$5 in a) print a[$1":"$2":"$4":"$5], exp($11)}' OFS="\t" "${G1000_EUR}.bim" "${INFILE}" > "${TEMP_APP1OUT}_pval.txt"

    #create sumstats file that is correctly formatted

    magma \
    --bfile $ref \
    --pval "${TEMP_APP1OUT}_pval.txt" N="${N}" \
    --gene-annot "${TEMP_APP2OUT}".annot \
    --out "${TEMP_APP3OUT}"

    #run pops
    python "${POPS}" \
    --gene_annot_path "${ANNOTFILE}" \
    --feature_mat_prefix "${MUNGE_FEATURES}" \
    --num_feature_chunks 2 \
    --magma_prefix "${TEMP_APP3OUT}" \
    --out_prefix "${OUTFILE_BASE}" \
    --verbose
}