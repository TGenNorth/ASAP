#!/usr/bin/env Rscript

library(tidyverse)
library(ASAPTools)

# Capture arguments passed from Nextflow
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: process_xml.R <xml_file> <min_snp> <sample_id>")
}

# Assign the arguments to variables
xml_file   <- args[1]
min_snp    <- as.numeric(args[2])
sample_id  <- args[3]

# 1. Individual Processing 
ASAP <- ASAPTools::read.ASAP.individual(xml_file)
SNPS <- ASAPTools::read.ASAP.snps.individual(xml_file)

# Define the columns that SHOULD be numeric
asap_numeric_names <- c("mapped_reads", "unassigned_reads", "unmapped_reads", "amplicon_number", "amplicon_reads", "breadth", "avg_depth")
snps_numeric_names <- c("Mapped_Reads", "unassigned_reads", "unmapped_reads", 
                        "amplicon_number", "location_depth", "snp_position", 
                        "snp_depth", "snp_proportion")

# 2. Numeric Cleaning (Using a safe helper function)
safe_convert <- function(df, target_cols) {
  # Only attempt to convert columns that actually exist in the data frame
  existing_cols <- intersect(target_cols, colnames(df))
  df[existing_cols] <- lapply(df[existing_cols], function(x) suppressWarnings(as.numeric(x)))
  return(df)
}

ASAP <- safe_convert(ASAP, asap_numeric_names)
SNPS <- safe_convert(SNPS, snps_numeric_names)

# 3. Filter SNPs (Removed the backslash \$ since this is a pure R file now)
SNPS <- SNPS[SNPS$snp_proportion >= min_snp, ]

cat("Line 38")

# 4. Extract Array Data
Depth <- ASAP.get.depth(ASAP, num_cores = 1)
Proportions <- ASAP.get.proportions(ASAP, num_cores = 1)
N_Reads <- ASAP.get.nreads(ASAP, num_cores = 1)
Quality.Discards <- ASAP.get.quality.discards(ASAP, num_cores = 1)

cat("Line 46")

# Joining results
array_info <- Depth %>%
  left_join(Proportions,     by = c("run", "name", "assay_name", "position")) %>% 
  left_join(N_Reads,         by = c("run", "name", "assay_name", "position")) %>%
  left_join(Quality.Discards, by = c("run", "name", "assay_name", "position"))

cat("Line 57")

# 5. Save outputs
save(ASAP, SNPS, array_info, file = paste0(sample_id, "_XML_Data.Rdata"))
write.csv(ASAP, file = paste0(sample_id, "_Summary.csv"), row.names = FALSE)