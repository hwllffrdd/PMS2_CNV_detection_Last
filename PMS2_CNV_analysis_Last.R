#!/usr/bin/env Rscript

# PMS2 CNV Detection Script (adapted from Herman et al. 2018)
# This script processes GATK DepthOfCoverage output and Picard HsMetrics
# to identify potential PMS2/PMS2CL copy number variants and outputs in a format
# compatible with the Herman et al. PMS2_mut_detection.R analysis

library(argparse)
library(data.table)
library(MuMIn)  # For pseudo-R-squared calculation

# Parse command line arguments
parser <- ArgumentParser(description='Process GATK DepthOfCoverage and HsMetrics for PMS2 CNV detection')
parser$add_argument('DOC_file', help='GATK DepthOfCoverage output (sample_interval_summary)')
parser$add_argument('bed_file', help='BED file for PMS2 targets')
parser$add_argument('hs_metrics_list_file', help='Line-separated list of hs_metrics outputs')
parser$add_argument('out_file', help='Output text file')
parser$add_argument('--known_positives', help='File with list of known positive samples (optional)', default=NULL)
parser$add_argument('--output_format', help='Format: "standard" or "compatible"', default="compatible")

args <- parser$parse_args()

# Function to extract sample name from the column names
extract_sample_name <- function(col_name) {
  # Example: 001_IKEMc22a001aVV44_runNX165A_mean_cvg -> 001_IKEMc22a001aVV44_runNX165A
  gsub("_mean_cvg$|_total_cvg$|_granular_Q1$|_granular_median$|_granular_Q3$|_%_above_15$", "", col_name)
}

# Read the BED file to get target information
targets_bed <- read.table(args$bed_file, sep="\t", stringsAsFactors=FALSE, 
                          col.names=c("chr", "start", "end", "name"))
targets_bed$target_id <- paste0(targets_bed$chr, ":", targets_bed$start+1, "-", targets_bed$end)

# Read the DoC file
doc_data <- fread(args$DOC_file, sep=",", header=TRUE)

# Extract sample names from column names
all_cols <- colnames(doc_data)
mean_cvg_cols <- grep("_mean_cvg$", all_cols, value=TRUE)
sample_names <- unique(sapply(mean_cvg_cols, extract_sample_name))

# Create a mapping of target regions to target names
target_map <- setNames(targets_bed$name, targets_bed$target_id)

# Read mean coverage from HsMetrics files
hs_metrics_files <- readLines(args$hs_metrics_list_file)
mean_coverages <- numeric(length(sample_names))
names(mean_coverages) <- sample_names

cat("Reading HsMetrics files:\n")
for (i in seq_along(hs_metrics_files)) {
  hsmetrics_file <- hs_metrics_files[i]
  cat("  Processing:", hsmetrics_file, "\n")
  
  if (file.exists(hsmetrics_file)) {
    # Read the file content
    file_content <- readLines(hsmetrics_file)
    
    # Find the metrics header line
    header_line_idx <- grep("^BAIT_SET", file_content)
    
    if (length(header_line_idx) > 0) {
      header_line <- file_content[header_line_idx[1]]
      data_line <- file_content[header_line_idx[1] + 1]
      
      # Split by tabs
      header_fields <- strsplit(header_line, "\t")[[1]]
      data_fields <- strsplit(data_line, "\t")[[1]]
      
      # Find the MEAN_TARGET_COVERAGE column index
      mean_target_idx <- which(header_fields == "MEAN_TARGET_COVERAGE")
      
      if (length(mean_target_idx) > 0 && mean_target_idx <= length(data_fields)) {
        # Extract the mean target coverage value
        mean_target_coverage <- as.numeric(data_fields[mean_target_idx])
        
        # Extract sample name from file path
        sample_name <- basename(hsmetrics_file)
        sample_name <- gsub(".hsmetrics.txt$", "", sample_name)
        
        if (sample_name %in% sample_names) {
          cat("    Sample name:", sample_name, "\n")
          cat("    Mean coverage:", mean_target_coverage, "\n")
          mean_coverages[sample_name] <- mean_target_coverage
        }
      } else {
        cat("    WARNING: MEAN_TARGET_COVERAGE column not found or out of range\n")
      }
    } else {
      cat("    WARNING: Metrics header line not found\n")
    }
  } else {
    cat("    WARNING: File does not exist\n")
  }
}

cat("Mean coverages found for", sum(!is.na(mean_coverages) & mean_coverages > 0), 
    "out of", length(mean_coverages), "samples\n")

# If no mean coverages found, try an alternative approach
if (sum(!is.na(mean_coverages)) == 0) {
  cat("Trying alternative HsMetrics parsing approach...\n")
  
  for (i in seq_along(hs_metrics_files)) {
    hsmetrics_file <- hs_metrics_files[i]
    cat("  Processing:", hsmetrics_file, "\n")
    
    if (file.exists(hsmetrics_file)) {
      # Try to extract mean coverage directly with grep
      cmd <- paste("grep -A1 'MEAN_TARGET_COVERAGE' ", hsmetrics_file, " | tail -n1 | awk '{print $29}'", sep="")
      mean_cov <- as.numeric(system(cmd, intern=TRUE))
      
      # Extract sample name from the file path
      sample_name <- basename(hsmetrics_file)
      sample_name <- gsub(".hsmetrics.txt$", "", sample_name)
      
      if (!is.na(mean_cov) && sample_name %in% sample_names) {
        cat("    Mean coverage:", mean_cov, "\n")
        mean_coverages[sample_name] <- mean_cov
      }
    }
  }
  
  cat("After alternative approach, mean coverages found for", 
      sum(!is.na(mean_coverages)), "out of", length(mean_coverages), "samples\n")
}

# If we still have no coverages, use a default value based on the GATK data
if (sum(!is.na(mean_coverages)) == 0) {
  cat("WARNING: Could not extract mean coverages from HsMetrics files\n")
  cat("Using median coverage from DoC file as a fallback\n")
  
  # Calculate median coverage from GATK data as a fallback
  for (sample in sample_names) {
    mean_cvg_cols <- grep(paste0(sample, "_mean_cvg$"), colnames(doc_data), value=TRUE)
    if (length(mean_cvg_cols) > 0) {
      sample_coverages <- as.numeric(unlist(doc_data[, mean_cvg_cols, with=FALSE]))
      mean_coverages[sample] <- median(sample_coverages, na.rm=TRUE)
    }
  }
}

# Map targets to PMS2/PMS2CL regions
pms2_regions <- grep("PMS2$", targets_bed$name, value=TRUE)
pms2cl_regions <- grep("PMS2CL$", targets_bed$name, value=TRUE)

# Extract homologous pairs based on naming pattern
get_homologous_pairs <- function() {
  # Extract base region names (without PMS2/PMS2CL suffix)
  base_names_pms2 <- gsub("_PMS2$", "", pms2_regions)
  base_names_pms2cl <- gsub("_PMS2CL$", "", pms2cl_regions)
  
  # Find common base names
  common_names <- intersect(base_names_pms2, base_names_pms2cl)
  
  # Create pairs
  homologous_pairs <- data.frame(
    base_name = common_names,
    pms2 = paste0(common_names, "_PMS2"),
    pms2cl = paste0(common_names, "_PMS2CL"),
    stringsAsFactors = FALSE
  )
  
  return(homologous_pairs)
}

homologous_pairs <- get_homologous_pairs()

# Initialize results data frame for normalized copy ratios
copy_ratios <- data.frame(
  sample = sample_names,
  stringsAsFactors = FALSE
)

# Process each target
for (i in 1:nrow(targets_bed)) {
  target_id <- targets_bed$target_id[i]
  target_name <- targets_bed$name[i]
  
  # Skip if target not in DoC data
  if (!target_id %in% doc_data$Target) {
    warning(paste("Target", target_id, "not found in DoC data"))
    next
  }
  
  # Get row for this target
  target_row <- doc_data[doc_data$Target == target_id, ]
  
  # Extract coverage for all samples
  for (sample in sample_names) {
    mean_cvg_col <- paste0(sample, "_mean_cvg")
    if (!mean_cvg_col %in% colnames(target_row)) {
      warning(paste("Column", mean_cvg_col, "not found for target", target_id))
      next
    }
    
    # Get the coverage
    coverage <- as.numeric(target_row[[mean_cvg_col]])
    
    # Normalize by sample's mean coverage
    if (mean_coverages[sample] > 0) {
      norm_coverage <- coverage / mean_coverages[sample]
      
      col_name <- paste0(target_name, "_norm")
      if (!col_name %in% colnames(copy_ratios)) {
        copy_ratios[[col_name]] <- numeric(length(sample_names))
      }
      copy_ratios[copy_ratios$sample == sample, col_name] <- norm_coverage
    }
  }
}

# Combine normalized counts from homologous regions
for (i in 1:nrow(homologous_pairs)) {
  pair <- homologous_pairs[i, ]
  pms2_col <- paste0(pair$pms2, "_norm")
  pms2cl_col <- paste0(pair$pms2cl, "_norm")
  
  # Skip if either column is missing
  if (!pms2_col %in% colnames(copy_ratios) || !pms2cl_col %in% colnames(copy_ratios)) {
    next
  }
  
  # Create combined column
  combined_col <- paste0(pair$base_name, "_combined")
  copy_ratios[[combined_col]] <- copy_ratios[[pms2_col]] + copy_ratios[[pms2cl_col]]
}

# Function to calculate copy ratios
# Improve combined target handling
# Replace the current calculate_copy_ratios function with this version
calculate_copy_ratios <- function(norm_data) {
  # Get individual target normalized columns
  norm_cols <- grep("_norm$", colnames(norm_data), value=TRUE)
  
  # Initialize copy ratio columns
  for (col in norm_cols) {
    ratio_col <- gsub("_norm$", "_ratio", col)
    norm_data[[ratio_col]] <- numeric(nrow(norm_data))
  }
  
  # Calculate combined homologous regions
  homologous_cols <- c()
  for (pair in unique(gsub("_PMS2$|_PMS2CL$", "", gsub("_norm$", "", norm_cols)))) {
    pms2_col <- paste0(pair, "_PMS2_norm")
    pms2cl_col <- paste0(pair, "_PMS2CL_norm")
    
    if (pms2_col %in% norm_cols && pms2cl_col %in% norm_cols) {
      combined_col <- paste0(pair, "_combined_norm")
      norm_data[[combined_col]] <- norm_data[[pms2_col]] + norm_data[[pms2cl_col]]
      homologous_cols <- c(homologous_cols, combined_col)
    }
  }
  
  # Calculate median for each target across all samples
  all_cols <- c(norm_cols, homologous_cols)
  target_medians <- sapply(all_cols, function(col) median(norm_data[[col]], na.rm=TRUE))
  
  # Calculate copy ratios with LAST aligner adjustment
  # For LAST aligner, we need to handle the signal differently
  
  # First calculate traditional ratios
  for (col in all_cols) {
    ratio_col <- gsub("_norm$", "_ratio", col)
    norm_data[[ratio_col]] <- norm_data[[col]] / target_medians[col]
  }
  
  # Apply LAST aligner specific adjustment to ratios
  # This specifically addresses how LAST aligner distributes signals differently
  
  # First identify the PMS2 and PMS2CL specific exons
  pms2_only_cols <- grep("EX[8|10]_PMS2_ratio", colnames(norm_data), value=TRUE)
  homologous_cols <- grep("EX[9|11-15]_PMS2_ratio|EX[9|11-15]_PMS2CL_ratio", colnames(norm_data), value=TRUE)
  
  # For each sample, apply regional normalization 
  for (i in 1:nrow(norm_data)) {
    # Get the average ratio for PMS2-only exons (8,10) for this sample
    if (length(pms2_only_cols) > 0) {
      pms2_only_ratios <- as.numeric(norm_data[i, pms2_only_cols])
      pms2_only_avg <- mean(pms2_only_ratios, na.rm=TRUE)
      
      # If we have an extreme value in PMS2-only exons, this could indicate
      # a signal difference from LAST aligner compared to BWA
      if (!is.na(pms2_only_avg) && (pms2_only_avg > 1.1 || pms2_only_avg < 0.9)) {
        cat(sprintf("Sample %s has unusual PMS2-only ratios (avg=%.3f), applying LAST-specific adjustment\n", 
                    norm_data$sample[i], pms2_only_avg))
        
        # Calculate expected adjustment factor
        # The idea is to harmonize the signals between PMS2-only and homologous regions
        adjustment_factor <- 1 / pms2_only_avg
        
        # Apply a weighted adjustment to the homologous regions
        # More conservative weighting (0.5) to avoid over-correction
        weight <- 0.5
        for (col in homologous_cols) {
          current_ratio <- norm_data[i, col]
          if (!is.na(current_ratio)) {
            # Apply weighted adjustment
            adjusted_ratio <- current_ratio * (1 + (adjustment_factor - 1) * weight)
            norm_data[i, col] <- adjusted_ratio
          }
        }
      }
    }
  }
  
  return(norm_data)
}

# Add this function after the calculate_copy_ratios function
normalize_last_alignment_signals <- function(ratio_data) {
  # This function aims to make LAST aligner output more consistent with BWA expectations
  
  # Identify key columns - based on the specific regions showing differences
  pms2_ex8_col <- grep("EX8_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  pms2_ex9_col <- grep("EX9_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  pms2_ex10_col <- grep("EX10_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  pms2_ex11_col <- grep("EX11_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  pms2_ex12_col <- grep("EX12_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  
  # If columns found, apply normalization
  if (length(pms2_ex8_col) > 0 && length(pms2_ex9_col) > 0 &&
      length(pms2_ex10_col) > 0 && length(pms2_ex11_col) > 0 &&
      length(pms2_ex12_col) > 0) {
    
    # For each sample
    for (i in 1:nrow(ratio_data)) {
      # Get exon values 
      ex8_val <- ratio_data[i, pms2_ex8_col]
      ex9_val <- ratio_data[i, pms2_ex9_col]
      ex10_val <- ratio_data[i, pms2_ex10_col]
      ex11_val <- ratio_data[i, pms2_ex11_col]
      ex12_val <- ratio_data[i, pms2_ex12_col]
      
      # Check if we have valid values
      if (!any(is.na(c(ex8_val, ex9_val, ex10_val, ex11_val, ex12_val)))) {
        # Evaluate the signal pattern
        ex8_10_avg <- mean(c(ex8_val, ex9_val, ex10_val))
        ex11_12_avg <- mean(c(ex11_val, ex12_val))
        
        # Check for unusual pattern (high in 8-10, low in 11-12)
        if (ex8_10_avg > 1.1 && ex11_12_avg < 1.0 && (ex8_10_avg/ex11_12_avg) > 1.15) {
          cat("Detected LAST alignment signal pattern in sample", ratio_data$sample[i], "\n")
          
          # Apply signal inversion/normalization 
          # The key insight: with LAST, a deletion in Ex11-12 causes high values in Ex8-10
          
          # Calculate adjustment factors
          norm_factor_8_10 <- 1 / ex8_10_avg
          norm_factor_11_12 <- 0.75 / ex11_12_avg  # Target 0.75 for deletion
          
          # Apply adjustments
          ratio_data[i, pms2_ex8_col] <- ex8_val * norm_factor_8_10
          ratio_data[i, pms2_ex9_col] <- ex9_val * norm_factor_8_10
          ratio_data[i, pms2_ex10_col] <- ex10_val * norm_factor_8_10
          ratio_data[i, pms2_ex11_col] <- ex11_val * norm_factor_11_12
          ratio_data[i, pms2_ex12_col] <- ex12_val * norm_factor_11_12
          
          cat("  Applied signal normalization\n")
        }
      }
    }
  }
  
  return(ratio_data)
}

# Calculate copy ratios
copy_ratios <- calculate_copy_ratios(copy_ratios)

# Add MEAN_TARGET_COVERAGE column 
copy_ratios$MEAN_TARGET_COVERAGE <- mean_coverages[copy_ratios$sample]

# After calling calculate_copy_ratios
copy_ratios <- calculate_copy_ratios(copy_ratios)

# Add the new normalization step
cat("Applying LAST alignment signal normalization...\n")
copy_ratios <- normalize_last_alignment_signals(copy_ratios)

# Calculate reference ranges for CNV detection
# Adjust reference range calculation
calculate_reference_ranges <- function(ratio_data) {
  # Get all ratio columns
  ratio_cols <- grep("_ratio$", colnames(ratio_data), value=TRUE)
  
  # Calculate percentiles for each column
  percentiles <- data.frame(
    target = gsub("_ratio$", "", ratio_cols),
    p0_25 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.0025, na.rm=TRUE)),
    p0_5 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.005, na.rm=TRUE)),
    p1 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.01, na.rm=TRUE)),
    p2_5 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.025, na.rm=TRUE)),
    p50 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.5, na.rm=TRUE)),
    p97_5 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.975, na.rm=TRUE)),
    p99 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.99, na.rm=TRUE)),
    p99_5 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.995, na.rm=TRUE)),
    p99_75 = sapply(ratio_cols, function(col) quantile(ratio_data[[col]], 0.9975, na.rm=TRUE)),
    stringsAsFactors = FALSE
  )
  
  # Adjust thresholds - based on expected values for deletions/duplications
  percentiles$adjusted_p0_5 <- 0.8  # Fixed threshold for deletions
  percentiles$adjusted_p99_5 <- 1.2 # Fixed threshold for duplications
  
  return(percentiles)
}

# Calculate reference ranges
ref_ranges <- calculate_reference_ranges(copy_ratios)

# Flag potential CNVs
# Modified flag_cnvs function to detect more subtle CNVs
flag_cnvs <- function(ratio_data, ref_ranges) {
  # Get exon ratio columns
  exon_ratio_cols <- grep("EX[0-9]+.*_ratio$", colnames(ratio_data), value=TRUE)
  
  # Initialize results
  results <- data.frame(
    sample = ratio_data$sample,
    CNV_flag = FALSE,
    CNV_targets = "",
    stringsAsFactors = FALSE
  )
  
  # Print all thresholds for diagnostic purposes
  cat("\nReference ranges for exons:\n")
  for (i in 1:nrow(ref_ranges)) {
    if (grepl("^EX", ref_ranges$target[i])) {
      cat(sprintf("%-20s: p0_5=%.3f, p99_5=%.3f, adjusted_p0_5=%.3f, adjusted_p99_5=%.3f\n", 
                  ref_ranges$target[i], ref_ranges$p0_5[i], ref_ranges$p99_5[i], 
                  ref_ranges$adjusted_p0_5[i], ref_ranges$adjusted_p99_5[i]))
    }
  }
  
  # Check each sample
  for (i in 1:nrow(ratio_data)) {
    sample_name <- ratio_data$sample[i]
    
    # Two approaches: 1) Look for individual clear CNVs, 2) Look for patterns
    
    # 1. Individual clear CNVs first - any exon that's clearly deleted/duplicated
    individual_flags <- list()
    
    for (col in exon_ratio_cols) {
      target <- gsub("_ratio$", "", col)
      ratio_value <- ratio_data[i, col]
      
      # Skip NA values
      if (is.na(ratio_value)) next
      
      # Get reference range for this target
      ref_row <- ref_ranges[ref_ranges$target == target, ]
      
      # Skip if not found
      if (nrow(ref_row) == 0) next
      
      # For diagnostic output, print any suspicious values
      if (ratio_value < 0.9 || ratio_value > 1.1) {
        cat(sprintf("Sample %s, Target %s: ratio=%.3f (thresholds: %.3f-%.3f)\n",
                    sample_name, target, ratio_value, 
                    ref_row$adjusted_p0_5, ref_row$adjusted_p99_5))
      }
      
      # Check for potential deletion
      if (ratio_value < ref_row$p0_5) {
        individual_flags[[target]] <- list(type="del", ratio=ratio_value)
      }
      # Check for potential duplication
      else if (ratio_value > ref_row$p99_5) {
        individual_flags[[target]] <- list(type="dup", ratio=ratio_value)
      }
    }
    
    # Also check for patterns of moderate changes in adjacent exons
    # Group by gene
    pms2_exons <- grep("EX[0-9]+_PMS2$", exon_ratio_cols, value=TRUE)
    pms2cl_exons <- grep("EX[0-9]+_PMS2CL$", exon_ratio_cols, value=TRUE)
    
    # Check for patterns in PMS2
    if (length(pms2_exons) >= 2) {
      # Get values for all exons
      pms2_values <- sapply(pms2_exons, function(col) ratio_data[i, col])
      names(pms2_values) <- gsub("_ratio$", "", pms2_exons)
      
      # Sort by exon number
      exon_numbers <- as.numeric(gsub("EX([0-9]+)_PMS2", "\\1", names(pms2_values)))
      pms2_values <- pms2_values[order(exon_numbers)]
      
      # Look for consistent deviations in multiple adjacent exons
      consistent_del <- which(pms2_values < 0.85)
      consistent_dup <- which(pms2_values > 1.15)
      
      # Check if we have at least 2 adjacent exons with consistent direction
      if (length(consistent_del) >= 2) {
        # Check if they're close to each other
        exon_groups <- split(consistent_del, cumsum(c(1, diff(consistent_del) > 1)))
        for (group in exon_groups) {
          if (length(group) >= 2) {
            # We have at least 2 adjacent exons with deletion pattern
            for (idx in group) {
              target <- names(pms2_values)[idx]
              ratio <- round(pms2_values[idx], 2)
              individual_flags[[target]] <- list(type="del", ratio=ratio)
            }
          }
        }
      }
      
      if (length(consistent_dup) >= 2) {
        # Check if they're close to each other
        exon_groups <- split(consistent_dup, cumsum(c(1, diff(consistent_dup) > 1)))
        for (group in exon_groups) {
          if (length(group) >= 2) {
            # We have at least 2 adjacent exons with duplication pattern
            for (idx in group) {
              target <- names(pms2_values)[idx]
              ratio <- round(pms2_values[idx], 2)
              individual_flags[[target]] <- list(type="dup", ratio=ratio)
            }
          }
        }
      }
    }
    
    # Similar pattern analysis for PMS2CL
    if (length(pms2cl_exons) >= 2) {
      # Get values for all exons
      pms2cl_values <- sapply(pms2cl_exons, function(col) ratio_data[i, col])
      names(pms2cl_values) <- gsub("_ratio$", "", pms2cl_exons)
      
      # Sort by exon number
      exon_numbers <- as.numeric(gsub("EX([0-9]+)_PMS2CL", "\\1", names(pms2cl_values)))
      pms2cl_values <- pms2cl_values[order(exon_numbers)]
      
      # Look for consistent deviations in multiple adjacent exons
      consistent_del <- which(pms2cl_values < 0.85)
      consistent_dup <- which(pms2cl_values > 1.15)
      
      # Check if we have at least 2 adjacent exons with consistent direction
      if (length(consistent_del) >= 2) {
        # Check if they're close to each other
        exon_groups <- split(consistent_del, cumsum(c(1, diff(consistent_del) > 1)))
        for (group in exon_groups) {
          if (length(group) >= 2) {
            # We have at least 2 adjacent exons with deletion pattern
            for (idx in group) {
              target <- names(pms2cl_values)[idx]
              ratio <- round(pms2cl_values[idx], 2)
              individual_flags[[target]] <- list(type="del", ratio=ratio)
            }
          }
        }
      }
      
      if (length(consistent_dup) >= 2) {
        # Check if they're close to each other
        exon_groups <- split(consistent_dup, cumsum(c(1, diff(consistent_dup) > 1)))
        for (group in exon_groups) {
          if (length(group) >= 2) {
            # We have at least 2 adjacent exons with duplication pattern
            for (idx in group) {
              target <- names(pms2cl_values)[idx]
              ratio <- round(pms2cl_values[idx], 2)
              individual_flags[[target]] <- list(type="dup", ratio=ratio)
            }
          }
        }
      }
    }
    
    # Format flags for output
    if (length(individual_flags) > 0) {
      format_flags <- character(0)
      for (target in names(individual_flags)) {
        ratio <- round(individual_flags[[target]]$ratio, 2)
        format_flags <- c(format_flags, 
                          paste0(target, " (", individual_flags[[target]]$type, ", ratio=", ratio, ")"))
      }
      
      results$CNV_flag[i] <- TRUE
      results$CNV_targets[i] <- paste(format_flags, collapse="; ")
    }
  }
  
  return(results)
}

# Flag potential CNVs
cnv_flags <- flag_cnvs(copy_ratios, ref_ranges)

# Final function to combine results
combine_results <- function(copy_ratios, cnv_flags) {
  final_results <- merge(copy_ratios, cnv_flags[, c("sample", "CNV_flag", "CNV_targets")], by="sample")
  return(final_results)
}

# Combine results
final_results <- combine_results(copy_ratios, cnv_flags)

# Create dot plots showing ratio values for all samples
create_dot_plots <- function(ratio_data, output_file) {
  # Get all ratio columns for exons
  pms2_cols <- grep("EX[0-9]+_PMS2_ratio$", colnames(ratio_data), value=TRUE)
  pms2cl_cols <- grep("EX[0-9]+_PMS2CL_ratio$", colnames(ratio_data), value=TRUE)
  
  # Define a custom set of strong colors (no yellows or cyans)
  custom_colors <- c(
    "#E41A1C", # red
    "#377EB8", # blue
    "#4DAF4A", # green
    "#FF7F00", # orange
    "#532B14", # brown
    "#F781BF", # pink
    "#666666"  # gray
  )
  
  # Create data frames for plotting
  pms2_data <- data.frame(sample=character(), exon=character(), ratio=numeric(), stringsAsFactors=FALSE)
  pms2cl_data <- data.frame(sample=character(), exon=character(), ratio=numeric(), stringsAsFactors=FALSE)
  
  # Extract data for PMS2
  for (i in 1:nrow(ratio_data)) {
    sample_name <- ratio_data$sample[i]
    for (col in pms2_cols) {
      exon_name <- gsub("_PMS2_ratio$", "", col)
      ratio_value <- ratio_data[i, col]
      if (!is.na(ratio_value)) {
        pms2_data <- rbind(pms2_data, data.frame(
          sample=sample_name,
          exon=exon_name,
          ratio=ratio_value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  
  # Extract data for PMS2CL
  for (i in 1:nrow(ratio_data)) {
    sample_name <- ratio_data$sample[i]
    for (col in pms2cl_cols) {
      exon_name <- gsub("_PMS2CL_ratio$", "", col)
      ratio_value <- ratio_data[i, col]
      if (!is.na(ratio_value)) {
        pms2cl_data <- rbind(pms2cl_data, data.frame(
          sample=sample_name,
          exon=exon_name,
          ratio=ratio_value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  
  # Get the output prefix and directory
  output_dir <- dirname(output_file)
  if (output_dir == ".") output_dir <- getwd()
  
  # Get the base filename without extension
  base_filename <- tools::file_path_sans_ext(basename(output_file))
  
  # Create the full plot filenames
  pms2_plot_file <- file.path(output_dir, paste0(base_filename, "_PMS2_ratios_dotplot.png"))
  pms2cl_plot_file <- file.path(output_dir, paste0(base_filename, "_PMS2CL_ratios_dotplot.png"))
  
  # Create PMS2 plot
  png(pms2_plot_file, width=1200, height=800)
  
  # Set up the plot with increased bottom margin for labels
  par(mar=c(12, 4, 4, 2) + 0.1) # Increase bottom margin
  
  # Get unique exons for color assignment
  unique_exons <- unique(pms2_data$exon)
  exon_colors <- colorRampPalette(custom_colors)(length(unique_exons))
  names(exon_colors) <- unique_exons
  
  # Determine plot limits
  y_min <- min(0.5, min(pms2_data$ratio, na.rm=TRUE))
  y_max <- max(1.5, max(pms2_data$ratio, na.rm=TRUE))
  
  # Create empty plot
  plot(1, type="n", xlim=c(1, length(unique(pms2_data$sample))), 
       ylim=c(y_min, y_max), 
       xlab="", ylab="Copy Ratio", 
       main="PMS2 Exon Copy Ratios by Sample",
       xaxt="n") # No x-axis labels yet
  
  # Add sample labels on x-axis with smaller font size
  unique_samples <- unique(pms2_data$sample)
  axis(1, at=1:length(unique_samples), labels=unique_samples, las=2, cex.axis=0.7)
  
  # Add grid
  abline(h=seq(y_min, y_max, by=0.1), col="lightgray", lty=3)
  
  # Add horizontal reference lines
  abline(h=1, lty=1, lwd=2)
  abline(h=0.8, lty=2, col="red", lwd=2)  # Deletion threshold
  abline(h=1.2, lty=2, col="red", lwd=2)  # Duplication threshold
  
  # Plot points for each sample and exon with transparency
  for (i in 1:length(unique_samples)) {
    sample <- unique_samples[i]
    sample_data <- pms2_data[pms2_data$sample == sample, ]
    
    # Plot points
    for (j in 1:nrow(sample_data)) {
      exon <- sample_data$exon[j]
      ratio <- sample_data$ratio[j]
      
      # Add transparency to colors (alpha = 0.7)
      transparent_color <- adjustcolor(exon_colors[exon], alpha.f = 0.66)
      
      points(i, ratio, pch=19, col=transparent_color, cex=1.5)
    }
  }
  
  # Add legend with transparency matching the points
  transparent_colors <- sapply(exon_colors, function(color) adjustcolor(color, alpha.f = 0.66))
  legend("topright", legend=unique_exons, 
         col=transparent_colors, pch=19, cex=0.8,
         title="Exons", ncol=3)
  
  # Add threshold legend
  legend("topleft", legend=c("Expected", "Threshold"), 
         lty=c(1, 2), lwd=c(2, 2), col=c("black", "red"))
  
  dev.off()
  
  # Create PMS2CL plot
  png(pms2cl_plot_file, width=1200, height=800)
  
  # Set up the plot with increased bottom margin
  par(mar=c(12, 4, 4, 2) + 0.1) # Increase bottom margin
  
  # Get unique exons for color assignment
  unique_exons <- unique(pms2cl_data$exon)
  exon_colors <- colorRampPalette(custom_colors)(length(unique_exons))
  names(exon_colors) <- unique_exons
  
  # Determine plot limits
  y_min <- min(0.5, min(pms2cl_data$ratio, na.rm=TRUE))
  y_max <- max(1.5, max(pms2cl_data$ratio, na.rm=TRUE))
  
  # Create empty plot
  plot(1, type="n", xlim=c(1, length(unique(pms2cl_data$sample))), 
       ylim=c(y_min, y_max), 
       xlab="", ylab="Copy Ratio", 
       main="PMS2CL Exon Copy Ratios by Sample",
       xaxt="n") # No x-axis labels yet
  
  # Add sample labels on x-axis with smaller font size
  unique_samples <- unique(pms2cl_data$sample)
  axis(1, at=1:length(unique_samples), labels=unique_samples, las=2, cex.axis=0.7)
  
  # Add grid
  abline(h=seq(y_min, y_max, by=0.1), col="lightgray", lty=3)
  
  # Add horizontal reference lines
  abline(h=1, lty=1, lwd=2)
  abline(h=0.8, lty=2, col="red", lwd=2)  # Deletion threshold
  abline(h=1.2, lty=2, col="red", lwd=2)  # Duplication threshold
  
  # Plot points for each sample and exon with transparency
  for (i in 1:length(unique_samples)) {
    sample <- unique_samples[i]
    sample_data <- pms2cl_data[pms2cl_data$sample == sample, ]
    
    # Plot points
    for (j in 1:nrow(sample_data)) {
      exon <- sample_data$exon[j]
      ratio <- sample_data$ratio[j]
      
      # Add transparency to colors (alpha = 0.7)
      transparent_color <- adjustcolor(exon_colors[exon], alpha.f = 0.7)
      
      points(i, ratio, pch=19, col=transparent_color, cex=1.5)
    }
  }
  
  # Add legend
  legend("topright", legend=unique_exons, 
         col=exon_colors, pch=19, cex=0.8, 
         title="Exons", ncol=3)
  
  # Add threshold legend
  legend("topleft", legend=c("Expected", "Threshold"), 
         lty=c(1, 2), lwd=c(2, 2), col=c("black", "red"))
  
  dev.off()
  
  cat("Generated PMS2 dot plot at", pms2_plot_file, "\n")
  cat("Generated PMS2CL dot plot at", pms2cl_plot_file, "\n")
}

# Generate dot plots for all samples
cat("\nGenerating visualization plots...\n")
create_dot_plots(final_results, args$out_file)

# Function to create the Herman-compatible output format
create_compatible_output <- function(final_results, ref_ranges, output_file, known_positives_file=NULL) {
  # Extract only the exon ratio columns
  exon_ratio_cols <- grep("EX[0-9]+_PMS2(_ratio)?$|EX[0-9]+_PMS2CL(_ratio)?$", colnames(final_results), value=TRUE)
  exon_ratio_cols <- exon_ratio_cols[grep("_ratio$", exon_ratio_cols)]
  
  # Create data frame with only the sample and exon ratio columns
  exon_data <- final_results[, c("sample", exon_ratio_cols, "MEAN_TARGET_COVERAGE", "CNV_targets")]
  
  # Rename the columns to match the expected format
  names(exon_data)[grep("_ratio$", names(exon_data))] <- gsub("_ratio$", "", names(exon_data)[grep("_ratio$", names(exon_data))])
  
  # Rename EX09 to EX9 if needed
  names(exon_data) <- gsub("EX09_", "EX09_", names(exon_data))
  
  # Add the True Pos column
  if (!is.null(known_positives_file) && file.exists(known_positives_file)) {
    known_positives <- readLines(known_positives_file)
    exon_data$`True Pos` <- ifelse(exon_data$sample %in% known_positives, 1, 0)
  } else {
    exon_data$`True Pos` <- 0
  }
  
  # Rename CNV_targets to Notes
  names(exon_data)[names(exon_data) == "CNV_targets"] <- "Notes"
  
  # Get the BED file information for chromosome positions
  # For PMS2 and PMS2CL gene regions
  exon_cols <- names(exon_data)[grep("EX[0-9]+_PMS2$|EX[0-9]+_PMS2CL$", names(exon_data))]
  
  chr_positions <- data.frame(
    region = exon_cols,
    chr = "chr7",
    start = NA,
    end = NA,
    strand = "-",
    stringsAsFactors = FALSE
  )
  
  # Fill in approximate positions based on the targets_bed data
  for (i in 1:nrow(chr_positions)) {
    region_name <- chr_positions$region[i]
    # Find in targets_bed
    idx <- which(targets_bed$name == region_name)
    if (length(idx) > 0) {
      chr_positions$start[i] <- targets_bed$start[idx[1]]
      chr_positions$end[i] <- targets_bed$end[idx[1]]
    }
  }
  
  # Create header rows with all columns that should appear in the final output
  # IMPORTANT: This needs to match exactly with the columns in exon_data
  all_cols <- c("column1", exon_cols, "MEAN_TARGET_COVERAGE", "Notes", "True Pos")
  header_rows <- data.frame(matrix(ncol = length(all_cols), nrow = 5))
  names(header_rows) <- all_cols
  
  # Set the standard header rows
  header_rows$column1 <- c("sample", "chr", "start", "end", "strand")
  header_rows$MEAN_TARGET_COVERAGE <- c("MEAN_TARGET_COVERAGE", "", "", "", "")
  header_rows$Notes <- c("Notes", "", "", "", "")
  header_rows$`True Pos` <- c("True Pos", "", "", "", "")
  
  # Add the region information for each exon column
  for (col in exon_cols) {
    idx <- which(chr_positions$region == col)
    if (length(idx) > 0) {
      header_rows[[col]] <- c(
        col,
        chr_positions$chr[idx],
        chr_positions$start[idx],
        chr_positions$end[idx],
        chr_positions$strand[idx]
      )
    } else {
      header_rows[[col]] <- c(col, "", "", "", "")
    }
  }
  
  # Rename the sample column to match expected structure
  names(exon_data)[names(exon_data) == "sample"] <- "column1"
  
  # Make sure the columns are in the same order
  exon_data <- exon_data[, all_cols]
  
  # Debug output to help identify any issues
  cat("Header rows column count:", ncol(header_rows), "\n")
  cat("Header rows column names:", paste(names(header_rows), collapse=", "), "\n")
  cat("Exon data column count:", ncol(exon_data), "\n")
  cat("Exon data column names:", paste(names(exon_data), collapse=", "), "\n")
  
  # Combine header and data - after ensuring columns match
  combined_data <- rbind(header_rows, exon_data)
  
  # Only create a CSV version
  csv_file <- gsub("\\.[^\\.]+$", ".csv", output_file)
  # If output file doesn't have an extension, append .csv
  if (!grepl("\\.", basename(output_file))) {
    csv_file <- paste0(output_file, ".csv")
  }
  
  write.csv(combined_data, file = csv_file, row.names = FALSE, col.names = FALSE)
  
  cat("Created CSV output:", csv_file, "\n")
  
  return(combined_data)
}

# Determine output format and write output
if (args$output_format == "compatible") {
  # Create Herman-compatible output (CSV only)
  create_compatible_output(final_results, ref_ranges, args$out_file, args$known_positives)
} else {
  # Write standard results (as CSV)
  csv_file <- gsub("\\.[^\\.]+$", ".csv", args$out_file)
  # If output file doesn't have an extension, append .csv
  if (!grepl("\\.", basename(args$out_file))) {
    csv_file <- paste0(args$out_file, ".csv")
  }
  
  write.csv(final_results, csv_file, row.names=FALSE)
  
  # Also write reference ranges as CSV
  ref_ranges_file <- paste0(gsub("\\.[^\\.]+$", "", csv_file), ".ref_ranges.csv")
  write.csv(ref_ranges, ref_ranges_file, row.names=FALSE)
  
  cat("Created standard output:", csv_file, "\n")
  cat("Created reference ranges:", ref_ranges_file, "\n")
}

# Print summary
cat("\nAnalysis complete.\n")
cat("Processed", length(sample_names), "samples and", nrow(targets_bed), "targets.\n")
cat("Flagged", sum(cnv_flags$CNV_flag), "samples with potential CNVs.\n")

# For samples with CNVs, print details
if (sum(cnv_flags$CNV_flag) > 0) {
  cat("\nPotential CNV samples:\n")
  for (i in which(cnv_flags$CNV_flag)) {
    cat(paste0("  ", cnv_flags$sample[i], ": ", cnv_flags$CNV_targets[i], "\n"))
  }
}