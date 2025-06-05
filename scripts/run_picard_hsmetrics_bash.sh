# Create a directory for results
mkdir -p ~/data/analysis/hsmetrics

# Create a file that will contain paths to all HsMetrics results
> ~/data/analysis/hsmetrics_list.txt

for bam in $(find ~/data/fastq/aligned_last -name "*.sorted.bam"); do
    sample=$(basename $(dirname $bam))
    
    echo "Processing ${sample}..."
    
    # Run Picard HsMetrics with the hybrid reference
    java -jar $PICARD_JAR CollectHsMetrics \
        I=$bam \
        O=~/data/analysis/hsmetrics/${sample}.hsmetrics.txt \
        R=$REF_GENOME \
        BAIT_INTERVALS=~/data/pms2_targets.interval_list \
        TARGET_INTERVALS=~/data/pms2_targets.interval_list \
        VALIDATION_STRINGENCY=LENIENT
    
    # Add path to list file
    echo "~/data/analysis/hsmetrics/${sample}.hsmetrics.txt" >> ~/data/analysis/hsmetrics_list.txt
done