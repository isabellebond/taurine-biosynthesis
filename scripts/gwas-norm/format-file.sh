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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_command() {
    # Decompose the line array, this is really for code clarity
    # Make sure order matches columns in jaf.txt file (ordered in selcol object in jaf_preparation.py)
    INFILE=($(echo "${PROCESS_LINE[1]}" | tr '|' ' '))
    OUTFILE="${PROCESS_LINE[2]}"

    info_msg "# infiles=""${#INFILE}"
    info_msg "idx (SGE_TASK_ID) [${IDX}] == rowidx [${ROW_IDX}]"
    
    # This generates a temp outfile name in out working directory that we
    # will then move to the final location upon successful completion of the
    # task
    OUTFILE_BASE="$(basename "$OUTFILE")"
    TEMP_OUTFILE="$WORKING_DIR"/"$OUTFILE_BASE"
    
    info_msg "output file basename=""$OUTFILE_BASE"
    info_msg "temporary output file=""$TEMP_OUTFILE"
    
#     # NOTE: CHECK if env set, otherwise use sensible default
#     # Allowing for unset variables
    set +u
    if [ -z "$PYTHON_INTERPRETER" ]; then
         PYTHON_INTERPRETER="${HOME}/miniconda3/envs/taurine/bin/python"
    fi
    set -u
    APPLICATION1="${HOME}/taurine-biosynthesis/scripts/gwas-norm/format-file.py"
    
    # If the input file doesn't exist, run the application
    "${PYTHON_INTERPRETER}" "${APPLICATION1}" --infile "${INFILE}" --outfile "${TEMP_OUTFILE}"
    info_msg ""${PYTHON_INTERPRETER}" "${APPLICATION1}" --infile "${INFILE}" --outfile "${TEMP_OUTFILE}""
    info_msg "Finished running application: ${APPLICATION1}"
    #Move tmp to outfile
    mv "${TEMP_OUTFILE}" "${OUTFILE}"

 }
 . task_setup.sh