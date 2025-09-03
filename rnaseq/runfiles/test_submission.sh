for CORES in {1..21..4}; do
    echo "Submitting job with $CORES cores..."
    qsub -pe smp $CORES run_benchmark_job.sh $CORES
done
