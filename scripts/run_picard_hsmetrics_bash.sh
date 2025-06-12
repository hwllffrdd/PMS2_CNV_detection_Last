WORKDIR="/home/rstudio/data"

# Create a directory for results
mkdir -p $WORKDIR/analysis/hsmetrics

# Create a file that will contain paths to all HsMetrics results
> $WORKDIR/analysis/hsmetrics_list.txt

for bam in $(find $WORKDIR/fastq/aligned_last -name "*.sorted.bam"); do
    sample=$(basename $(dirname $bam))
    
    echo "Processing ${sample}..."
    
    # Run Picard HsMetrics with the hybrid reference
    java -jar $PICARD_JAR CollectHsMetrics \
        I=$bam \
        O=$WORKDIR/analysis/hsmetrics/${sample}.hsmetrics.txt \
        R=$REF_GENOME \
        BAIT_INTERVALS=$WORKDIR/pms2_targets.interval_list \
        TARGET_INTERVALS=$WORKDIR/pms2_targets.interval_list \
        VALIDATION_STRINGENCY=LENIENT
    
    # Add path to list file
    echo "$WORKDIR/analysis/hsmetrics/${sample}.hsmetrics.txt" >> $WORKDIR/analysis/hsmetrics_list.txt
done
