# Create a directory for results
mkdir -p ~/data/analysis

# Create a list of all BAM files
find ~/data/fastq/aligned_last -name "*.sorted.bam" > ~/data/analysis/all_bams.list

# Run GATK DepthOfCoverage on all samples
gatk --java-options "-Xmx8g" DepthOfCoverage \
    -R $REF_GENOME \
    -I ~/data/analysis/all_bams.list \
    -L ~/data/pms2_targets.bed \
    -O ~/data/analysis/combined_coverage \
    --omit-depth-output-at-each-base \
    --include-ref-n-sites