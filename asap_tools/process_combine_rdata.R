#!/usr/bin/env Rscript

# Load necessary libraries
# tidyverse for data manipulation, data.table for high-speed binding
library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
library(parallelly)

# 1. Capture Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: process_combine_rdata.R <optional_poi_csv> <input_files...>", call. = FALSE)
}

poi_csv <- args[1] # Will be "NULL" if not provided
files <- args[2:length(args)]

# Test env
# files <- list.files("/scratch/tporter/ASAP_TB_Validation/ASAP_TB_Subset/XML_Rdata/", pattern = "Rdata", full.names = T)
# poi_csv <- "/scratch/tporter/ASAP_TB_Validation/Updated_TB_Genes.csv"

# 2. Setup Parallel Backend
# parallelly::availableCores() is SLURM-aware and respects cpus allocated to the job
num_cores <- parallelly::availableCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

message(paste("🚀 Combining", length(files), "samples using", num_cores, "cores..."))

# 3. Parallel Loading with foreach
# We use %dopar% to read files simultaneously across available cores
combined_list <- foreach(f = files, .packages = c("tidyverse")) %dopar% {
  
  if (!file.exists(f)) {
    stop(paste("File not found in task directory:", f))
  }
  
  # Create a local environment for each file to prevent object collision
  temp_env <- new.env()
  load(f, envir = temp_env)
  
  # 2. Handle Positions of Interest (Optional)
  if (is.na(poi_csv) || poi_csv == "NULL" || poi_csv == "") {
    message("No Positions of Interest provided. Generating reference data...")

    } else {
    message(paste("Loading positions of interest from:", poi_csv))
    genes <- read.csv(poi_csv)
    
    # Generate positions for each gene in the CSV
    Gene_Positions <- genes %>%
      rowwise() %>%
      do(data.frame(
        position = seq(min(.$start, .$end), max(.$start, .$end)), 
        gene = .$gene,
        reference = .$seqnames
      )) %>%
      ungroup()
    
    # filter unneeded array info...
    temp_env$array_info <- temp_env$array_info %>%
      semi_join(Gene_Positions, by = c("position" = "position", "assay_name" = "reference"))
    
    # Remove array_info fields from ASAP
    temp_env$ASAP$n_reads <- NULL
    temp_env$ASAP$quality_discards <- NULL
    temp_env$ASAP$proportions <- NULL
    temp_env$ASAP$depths <- NULL
    
    # Remove unneeded data from fields.
    # For TB Example, no way to easily clean SNPs to match POI filtering... 
  }
  
  
  # Return as a structured list for easier extraction
  list(
    asap = temp_env$ASAP,
    snps = temp_env$SNPS,
    info = temp_env$array_info
  )
}

# Explicitly stop the cluster to free system resources
stopCluster(cl)

# 4. Fast Binding with data.table
# rbindlist is written in C and is significantly faster than map_df or rbind
message("📊 Merging data frames...")

# Merge ASAP summary data
final_asap  <- as.data.frame(data.table::rbindlist(map(combined_list, "asap"), fill = TRUE))
gc() # Clear memory immediately after large merges

# Merge SNP data
final_snps  <- as.data.frame(data.table::rbindlist(map(combined_list, "snps"), fill = TRUE))
gc()

# Merge large Array Info (Depth/Proportions)
final_array <- as.data.frame(data.table::rbindlist(map(combined_list, "info"), fill = TRUE))
gc()

# 5. Save Combined Outputs
# Saving both as a compressed Rdata object and a flat CSV summary
message(paste("💾 Saving results"))

save(final_asap, final_snps, final_array, file = "Combined_ASAP_Data.Rdata")
write.csv(final_asap, file = "Combined_Summary.csv", row.names = FALSE)

message("✅ Success: Combined data saved to current working directory.")