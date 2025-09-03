#!/bin/bash -l

# Request one week of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=168:00:00

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=4G

# Set the name of the job.
#$ -N rnaseq

#$ -cwd
#$ -V
#$ -o /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/jobs/rnaseq.o
#$ -e /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/jobs/rnaseq.e
#$ -pe smp 1
module load java/temurin-17/17.0.2_8
export NXF_OPTS='-Xms1g -Xmx3g'
nextflow run nf-core/rnaseq \
    --input /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/input/samplesheet.csv \
    --outdir /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/taurine \
    --email rmgpibo@ucl.ac.uk \
    --gtf /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/reference/Homo_sapiens.GRCh38.114.gtf.gz \
    --fasta /home/rmgpibo/Scratch/taurine-biosynthesis/data/rnaseq/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
    -profile ucl_myriad \
    -resume
    