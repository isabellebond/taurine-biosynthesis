for CORES in {1..21..4}; do
    echo "Submitting job with $CORES cores..."
    qsub runfiles/test_cores.sh $CORES
done
