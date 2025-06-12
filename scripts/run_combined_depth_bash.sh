WORKDIR="/home/rstudio/data"
mkdir -p $WORKDIR/analysis
find $WORKDIR/fastq/aligned_last -name "*.sorted.bam" > $WORKDIR/analysis/all_bams.list
gatk --java-options "-Xmx8g" DepthOfCoverage \
    -R $REF_GENOME \
    -I $WORKDIR/analysis/all_bams.list \
    -L $WORKDIR/pms2_targets.bed \
    -O $WORKDIR/analysis/combined_coverage \
    --omit-depth-output-at-each-base \
    --include-ref-n-sites
