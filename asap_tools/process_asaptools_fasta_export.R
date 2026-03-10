#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Check usage...")
}

rdata_input <- args[1]
prefix      <- args[2]
BREADTH_COVERAGE_THRESHOLD <- as.numeric(args[3])*100

# rdata_input <- "/scratch/tporter/ASAP_SC2_Validation/ASAP_Illumina_Paired_SE_Test/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# poi_csv     <- "NULL"
#rdata_input <- "/scratch/tporter/ASAP_TB_Validation/ASAP_TB_Subset/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "/scratch/tporter/ASAP_TB_Validation/Updated_TB_Genes.csv"

unique_assays <- unique(final_asap$assay_name)

for(ASSAY in unique_assays) {

  # Filter data for the specific assay and threshold
  Fasta_df <- final_asap %>%
    filter(assay_name == ASSAY) %>%
    filter(breadth > BREADTH_COVERAGE_THRESHOLD)

  # Check if Fasta_df is empty
  if (nrow(Fasta_df) == 0) {
    message(paste("Skipping", ASSAY, "- No samples met the breadth threshold of", BREADTH_COVERAGE_THRESHOLD))
    next # Move to the next ASSAY in the loop
  }

  # Vectorized creation of FASTA lines: much faster than a loop
  # Format: >SampleName_AssayName \n Sequence
  fasta_lines <- paste0(">", Fasta_df$name, "_", ASSAY, "\n", Fasta_df$consensus_seq)

  # Clean the assay name for a safe filename
  clean_assay <- gsub("[^[:alnum:]]", "_", ASSAY)
  file_name <- paste0("./", prefix, "_", clean_assay, ".fasta")

  # Write the file
  writeLines(fasta_lines, file_name)

  message(paste("Successfully exported", nrow(Fasta_df), "sequences to", file_name))
}
