#!/bin/bash
#$ -l mem=8G             # Request 8G of memory *per core*
#$ -l h_rt=24:00:00      # Request 48 hours of wall time, matching the max
#$ -cwd                  # Run the job from the current working directory
#$ -j y                  # Merge standard error and standard out
#$ -S /bin/bash          # Specify the shell
#$ -N rnaseq_benchmark   # Job name
#$ -o /home/rmgpibo/Scratch/taurine-biosynthesis/rnaseq/jobs/rnaseq.o
#$ -e /home/rmgpibo/Scratch/taurine-biosynthesis/rnaseq/jobs/rnaseq.e

# --- Safety Check ---
# This ensures the script exits if a number of cores isn't provided
if [ -z "$1" ]; then
    echo "Error: Please provide the number of cores as the first argument."
    echo "Usage: qsub -pe smp <cores> run_benchmark_job.sh <cores>"
    exit 1
fi

CORES=$1

echo "=========================================================="
echo "Starting nf-core/rnaseq benchmark..."
echo "Job ID: $JOB_ID"
echo "Running on host: $(hostname)"
echo "Requested Cores: $CORES"
echo "=========================================================="

# --- Environment Setup ---
module load java/temurin-17/17.0.2_8
export NXF_OPTS='-Xms1g -Xmx3g'

# --- Define File Paths ---
# Using absolute paths is a good practice
BASE_DIR="/home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq"
FASTQ_DIR="${BASE_DIR}/input"
OUTDIR="${BASE_DIR}/benchmark"
REFERENCE_DIR="${BASE_DIR}/reference"
RUNFILES_DIR="${BASE_DIR}/runfiles"

# --- Run Nextflow Pipeline ---
nextflow run nf-core/rnaseq \
    -profile ucl_myriad \
    --input "${FASTQ_DIR}/samplesheet.csv" \
    --outdir "${OUTDIR}/test_${CORES}cores" \
    --gtf "${REFERENCE_DIR}/Homo_sapiens.GRCh38.114.gtf.gz" \
    --fasta "${REFERENCE_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" \
    -c "${RUNFILES_DIR}/benchmark.cnf" \
    --skip_qc \

echo "=========================================================="
echo "Benchmark with $CORES cores finished."
echo "=========================================================="
