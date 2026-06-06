#!/usr/bin/env Rscript

library(tidyverse)
library(xml2)

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_read.ASAP.individual.R"))
source(file.path(.functions_dir, "_read.ASAP.snps.individual.R"))
source(file.path(.functions_dir, "_ASAP.get.depth.R"))
source(file.path(.functions_dir, "_ASAP.get.proportions.R"))
source(file.path(.functions_dir, "_ASAP.get.nreads.R"))
source(file.path(.functions_dir, "_ASAP.get.quality.discards.R"))

# Capture arguments passed from Nextflow
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: process_xml.R <xml_file> <min_snp> <sample_id>")
}

# Assign the arguments to variables
xml_file   <- args[1]
min_snp    <- as.numeric(args[2])*100
sample_id  <- args[3]

# xml_file   <- "/scratch/tporter/ASAP_SC2_Results/ASAP_Illumina_Reads_R_Out_Sunday/xml/DRR640392.xml"
# min_snp    <- 0.01*100
# sample_id  <- "DRR640392"


# 1. Individual Processing
ASAP <- read.ASAP.individual(xml_file)
SNPS <- read.ASAP.snps.individual(xml_file)

# Define the columns that SHOULD be numeric
asap_numeric_names <- c("total_reads", "trimmed_reads", "mapped_reads", "unassigned_reads", "unmapped_reads",
                        "amplicon_number", "amplicon_reads", "aligned_reads",
                        "primer_reads", "no_primer_reads",
                        "identity_input", "identity_discarded",
                        "smor_input", "smor_pairs_dropped", "smor_consensus_reads",
                        "breadth", "avg_depth")
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
