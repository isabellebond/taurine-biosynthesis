#!/bin/bash
#$ -l mem=8G
#$ -l h_rt=24:00:00
#$ -pe smp 16
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N rnaseq_benchmark

module load java/temurin-17/17.0.2_8
export NXF_OPTS='-Xms1g -Xmx3g'

FASTQ_DIR="/home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/input"
OUTDIR="data/rnaseq/benchmark"
GENOME="GRCh38"

for CORES in 1 4 8 16; do
    echo "=== Running nf-core/rnaseq with $CORES cores ==="

    nextflow run nf-core/rnaseq \
        -profile singularity \
        --input "${FASTQ_DIR}/samplesheet.csv" \
        --outdir "${OUTDIR}/test_${CORES}cores" \
        --gtf /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/reference/Homo_sapiens.GRCh38.114.gtf.gz \
        --fasta /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
        -c /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/runfiles/benchmark.cnf \
        -profile ucl_myriad \
        -process.cpus $CORES \
        -with-report report_${CORES}cores.html \
        -with-trace trace_${CORES}cores.txt \
        -with-timeline timeline_${CORES}cores.html \
        -with-dag flowchart_${CORES}cores.png \
        -resume \
        --skip_qc \
        --max_reads 5000000 \
        -process.Name:STAR_ALIGN.cpus $CORES
done
