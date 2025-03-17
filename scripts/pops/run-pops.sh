#!/bin/bash
# Source in the task setup, this will deal with the command line input
# and cleanup functions, there are several variables that it will give
# access to:
# JA="$1"
# TMP_LOC="$2"
# STEP=$3
# IDX=${4:-$SGE_TASK_ID}

# WORKING_DIR : ~/Scratch/tmp/

# Store all the paths to the intermediate files, these will all be deleted when
# exitting cleanly or with an error. The "" is to stop unbound errors
# TEMPFILES=("")

# Files that should be moved should go into these two arrays. This will then be
# moved to their final location on exit.
# MOVE_FROM=()
# MOVE_TO=()

# PROCESS_LINE : An array extracted from the job array file, this provides
# input for the current STEP of the starting from the current TASK

# ROW_IDX : Element 0 from the process line, this is a requirement that all
# job arrays have this


# Enable nounset option
set -u
run_command() {
    # Decompose the line array, this is really for code clarity
    # Make sure order matches columns in jaf.txt file (ordered in selcol object in jaf_preparation.py)
    INFILE="${PROCESS_LINE[1]}" #GWAS sumstats in line 1
    N="${PROCESS_LINE[2]}" #Number of samples in line 2
    ANNOTFILE="${PROCESS_LINE[3]}"
    PFILE="${PROCESS_LINE[4]}"
    OUTFILE="${PROCESS_LINE[5]}" #Output file for analysis

    set +u
    if [ -z "$PYTHON_INTERPRETER" ]; then
         PYTHON_INTERPRETER="${HOME}/miniconda3/envs/taurine/bin/python"
    fi
    set -u

    MUNGE_FEATURES="/lustre/projects/mol_cardio/resources/software/pops/example/data/features_munged"
    GENE_ANNOT="/lustre/projects/mol_cardio/resources/software/pops/example/data/utils/gene_annot_jun10.txt"
    MAKE_SNPLOC="/home/rmgpibo/taurine-biosynthesis/scripts/pops/make_snploc.py"
    POPS="/lustre/projects/mol_cardio/resources/software/pops/pops.py"
    NCBI38="/lustre/projects/mol_cardio/resources/data/magma/NCBI38.gene.loc"
    G1000_EUR="/lustre/projects/mol_cardio/resources/data/magma/g1000_eur"

    

    # Create SNP location file using awk
    python "${MAKE_SNPLOC}" \
            --infile "${INFILE}" \
            --bimfile "${G1000_EUR}" \
            --outfile "${PFILE}_sumstats.txt.gz"

    if [ ! -f "${ANNOTFILE}.genes.annot" ]; then
        zcat "${ANNOTFILE}"_sumstats.txt.gz | awk 'NR>1 {print $1,$3,$4}' > "${ANNOTFILE}_snploc.txt"

        magma --annotate \
            --snp-loc "${ANNOTFILE}_snploc.txt" \
            --gene-loc "${NCBI38}" \
            --out "${ANNOTFILE}"
    else
        echo "Annotation file already exists: ${ANNOTFILE}.genes.annot"
    fi

    #Get snp and p value pairs for analysis
    zcat "${ANNOTFILE}"_sumstats.txt | awk 'NR>1 {print $1,$2}' > "${PFILE}_pvals.txt"
    #create sumstats file that is correctly formatted

    magma \
    --bfile $ref \
    --pval "${PFILE}_pvals.txt" N=14759 \
    --gene-annot "${ANNOTFILE}.annot" \
    --out "${OUTFILE_PREFIX}"

    #run pops
    python "${POPS}" \
    --gene_annot_path "${ANNOTFILE}" \
    --feature_mat_prefix "${MUNGE_FEATURES}" \
    --num_feature_chunks 2 \
    --magma_prefix "${TEMP_APP3OUT}" \
    --out_prefix "${OUTFILE_BASE}" \
    --verbose
}
. task_setup.sh