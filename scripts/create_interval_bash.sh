# Create header for interval list
samtools view -H ~/data/fastq/aligned_last/$(ls ~/data/fastq/aligned_last | head -1)/$(ls ~/data/fastq/aligned_last | head -1).sorted.bam | grep '@SQ' > header.txt
echo '@HD	VN:1.6	SO:coordinate' | cat - header.txt > interval_list_header.txt

# Convert BED to interval list - note the fixed format here
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$3,"+","TARGET_" NR}' ~/data/pms2_targets.bed > intervals.txt
cat interval_list_header.txt intervals.txt > ~/data/pms2_targets.interval_list
rm header.txt interval_list_header.txt intervals.txt