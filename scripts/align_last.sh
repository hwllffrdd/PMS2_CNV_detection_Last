#!/bin/bash
# Optimized LAST alignment script for PMS2/PMS2CL analysis

reference=$REF_GENOME
out_dir="aligned_last"
db_dir="${out_dir}/reference_db"

# Create main output and database directories
mkdir -p "$out_dir"
mkdir -p "$db_dir"

# Create LAST database once at the beginning
echo "Creating LAST database (once for all samples)..."
# Copy the reference FASTA to a writable location
cp "$reference" "${db_dir}/reference.fa"
samtools faidx "${db_dir}/reference.fa"
# Create LAST database with parameters for highly homologous regions
lastdb -uNEAR -R01 "${db_dir}/ref_db" "${db_dir}/reference.fa"

# Create a custom header directly using the reference
echo "Creating common header template..."
echo -e "@HD\tVN:1.6\tSO:coordinate" > "${db_dir}/header_template.sam"
# Extract sequence info from the reference to create SQ lines
cut -f1,2 "${db_dir}/reference.fa.fai" | awk '{print "@SQ\tSN:"$1"\tLN:"$2}' >> "${db_dir}/header_template.sam"

process_sample() {
    r1_file=$1
    sample_name=$(basename $r1_file _R1.fastq.gz)
    r2_file="${r1_file/_R1/_R2}"
    
    echo "Processing sample: $sample_name"
    mkdir -p "${out_dir}/${sample_name}"
    work_dir="${out_dir}/${sample_name}"
    
    # Align reads with LAST
    echo "  Aligning with LAST..."
    # For read 1
    lastal -Q1 -r1 -a15 -b3 "${db_dir}/ref_db" "$r1_file" > "${work_dir}/${sample_name}_r1.maf"
    # For read 2
    lastal -Q1 -r1 -a15 -b3 "${db_dir}/ref_db" "$r2_file" > "${work_dir}/${sample_name}_r2.maf"
    
    # Create sample-specific header by adding read group info to the template
    echo "  Creating sample-specific header..."
    cp "${db_dir}/header_template.sam" "${work_dir}/header.sam"
    echo -e "@RG\tID:$sample_name\tSM:$sample_name\tLB:lib1\tPL:illumina" >> "${work_dir}/header.sam"
    
    # Convert MAF to SAM with proper formatting
    echo "  Converting MAF to SAM..."
    maf-convert sam "${work_dir}/${sample_name}_r1.maf" | grep -v "^@" > "${work_dir}/${sample_name}_r1_body.sam"
    maf-convert sam "${work_dir}/${sample_name}_r2.maf" | grep -v "^@" > "${work_dir}/${sample_name}_r2_body.sam"
    
    # Add RG tag to each alignment line
    echo "  Adding read group tags to alignments..."
    awk -v rg="$sample_name" '{print $0"\tRG:Z:"rg}' "${work_dir}/${sample_name}_r1_body.sam" > "${work_dir}/${sample_name}_r1_rg_body.sam"
    awk -v rg="$sample_name" '{print $0"\tRG:Z:"rg}' "${work_dir}/${sample_name}_r2_body.sam" > "${work_dir}/${sample_name}_r2_rg_body.sam"
    
    # Combine header and aligned reads
    cat "${work_dir}/header.sam" "${work_dir}/${sample_name}_r1_rg_body.sam" > "${work_dir}/${sample_name}_r1.sam"
    cat "${work_dir}/header.sam" "${work_dir}/${sample_name}_r2_rg_body.sam" > "${work_dir}/${sample_name}_r2.sam"
    
    # Convert to BAM directly
    echo "  Converting to BAM format..."
    samtools view -b "${work_dir}/${sample_name}_r1.sam" > "${work_dir}/${sample_name}_r1.bam"
    samtools view -b "${work_dir}/${sample_name}_r2.sam" > "${work_dir}/${sample_name}_r2.bam"
    
    # Merge BAM files
    echo "  Merging BAM files..."
    samtools merge -f "${work_dir}/${sample_name}.bam" \
        "${work_dir}/${sample_name}_r1.bam" \
        "${work_dir}/${sample_name}_r2.bam"
    
    # Sort the merged BAM
    echo "  Sorting the merged BAM file..."
    samtools sort -o "${work_dir}/${sample_name}.sorted.bam" "${work_dir}/${sample_name}.bam"
    
    # Index BAM file
    echo "  Indexing the BAM file..."
    samtools index "${work_dir}/${sample_name}.sorted.bam"
    
    # Clean up intermediate files
    echo "  Cleaning up intermediate files..."
    rm "${work_dir}/${sample_name}_r1.maf" "${work_dir}/${sample_name}_r2.maf" \
       "${work_dir}/${sample_name}_r1_body.sam" "${work_dir}/${sample_name}_r2_body.sam" \
       "${work_dir}/${sample_name}_r1_rg_body.sam" "${work_dir}/${sample_name}_r2_rg_body.sam" \
       "${work_dir}/${sample_name}_r1.sam" "${work_dir}/${sample_name}_r2.sam" \
       "${work_dir}/${sample_name}_r1.bam" "${work_dir}/${sample_name}_r2.bam" \
       "${work_dir}/${sample_name}.bam" "${work_dir}/header.sam"
    
    echo "  Done processing ${sample_name}"
}

export -f process_sample
export reference
export out_dir
export db_dir

# Process all samples
if [ "$#" -eq 0 ]; then
    # If no arguments provided, process all R1 files in current directory
    for r1_file in *_R1.fastq.gz; do
        process_sample "$r1_file"
    done
else
    # Process specific samples provided as arguments
    for r1_file in "$@"; do
        process_sample "$r1_file"
    done
fi

echo "Cleaning up reference database..."
rm -rf "$db_dir"
echo "Database cleanup complete."
