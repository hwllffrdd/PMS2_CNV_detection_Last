# PMS2/PMS2CL Copy Number Variant Detection Pipeline

A bioinformatics pipeline for detecting copy number variants (CNVs) in the highly homologous PMS2/PMS2CL genes using LAST aligner and read depth analysis.

## Overview

This pipeline implements the methodology described in Herman et al. 2018 for detecting deletions and duplications in PMS2 exons 9-15, which are highly homologous to the PMS2CL pseudogene.

## Requirements

- Docker
- FASTQ files (paired-end, named with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes)
- At least 24GB RAM recommended

## Quick Start

1. **Build the Docker image:**
   ```bash
   docker build -f DockerfileLAST -t pms2_pipeline .
   ```

2. **Run the container:**
   ```bash
   docker run -v $(pwd):/home/rstudio/data --rm -it -m 24g --memory-swap 24g pms2_pipeline bash
   ```
   or with the use of RStudio Sever:
   ```bash
   docker run -v $(pwd):/home/rstudio/data --rm -e DISABLE_AUTH=true -m 24g --memory-swap 24g -p 8787:8787 pms2_pipeline
   ```
3. **Navigate to your FASTQ directory and run the pipeline:**
   ```bash
   cd /home/rstudio/data/fastq
   ```

## Pipeline Steps

### Step 1: Alignment
```bash
../scripts/align_last.sh
```
- Aligns FASTQ files to hg19 chromosome 7 using LAST aligner
- Creates sorted and indexed BAM files in `aligned_last/` directory
- Optimized for highly homologous sequences

### Step 2: Calculate Read Depth
```bash
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
```

### Step 3: Create Interval List
```bash
samtools view -H $WORKDIR/fastq/aligned_last/$(ls $WORKDIR/fastq/aligned_last | head -1)/$(ls $WORKDIR/fastq/aligned_last | head -1).sorted.bam | grep '@SQ' > header.txt
echo '@HD	VN:1.6	SO:coordinate' | cat - header.txt > interval_list_header.txt
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$3,"+","TARGET_" NR}' $WORKDIR/pms2_targets.bed > intervals.txt
cat interval_list_header.txt intervals.txt > $WORKDIR/pms2_targets.interval_list
rm header.txt interval_list_header.txt intervals.txt
```

### Step 4: Run Picard HsMetrics
```bash
mkdir -p $WORKDIR/analysis/hsmetrics
> $WORKDIR/analysis/hsmetrics_list.txt

for bam in $(find $WORKDIR/fastq/aligned_last -name "*.sorted.bam"); do
    sample=$(basename $(dirname $bam))
    echo "Processing ${sample}..."
    java -jar $PICARD_JAR CollectHsMetrics \
        I=$bam \
        O=$WORKDIR/analysis/hsmetrics/${sample}.hsmetrics.txt \
        R=$REF_GENOME \
        BAIT_INTERVALS=$WORKDIR/pms2_targets.interval_list \
        TARGET_INTERVALS=$WORKDIR/pms2_targets.interval_list \
        VALIDATION_STRINGENCY=LENIENT
    echo "$WORKDIR/analysis/hsmetrics/${sample}.hsmetrics.txt" >> $WORKDIR/analysis/hsmetrics_list.txt
done
```

### Step 5: CNV Analysis
```bash
mkdir -p $WORKDIR/formatted_results

Rscript PMS2_CNV_analysis_Last.R \
    $WORKDIR/analysis/combined_coverage.sample_interval_summary \
    $WORKDIR/pms2_targets.bed \
    $WORKDIR/analysis/hsmetrics_list.txt \
    $WORKDIR/formatted_results/result.csv \
    --known_positives $WORKDIR/known_positives.txt
# note: --known_positives is optional
```

## Output

Results are generated in the `formatted_results/` directory:
- CNV calls for each sample
- Copy ratio plots
- Summary statistics
- Flagged samples requiring follow-up

## Target Regions

The pipeline analyzes 38 target regions:
- 8 regions in PMS2 (exons 8-15)
- 6 homologous regions in PMS2CL
- Based on coordinates from Herman et al. 2018 (Table 1)

## References

Herman DS, Smith C, Liu C, et al. Efficient Detection of Copy Number Mutations in PMS2 Exons with a Close Homolog. J Mol Diagn. 2018;20(4):512-521.

## Files Included

- `DockerfileLAST` - Docker container definition
- `pms2_targets.bed` - Target regions for analysis
- `scripts/align_last.sh` - LAST alignment script
- `scripts/run_combined_depth_bash.sh` - GATK DepthOfCoverage commands
- `scripts/create_interval_bash.sh` - Interval list creation commands
- `scripts/run_picard_hsmetrics_bash.sh` - Picard HsMetrics commands
- `PMS2_CNV_analysis_Last.R` - R script for CNV analysis
- `known_positives.txt` - List of known positive control samples
