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

# 2. Numeric Cleaning 
cols_asap <- c(3, 4, 5, 11, 13, 14)
cols_snps <- c(8, 10, 12, 13)
ASAP[cols_asap] <- lapply(ASAP[cols_asap], function(x) suppressWarnings(as.numeric(x)))
SNPS[cols_snps] <- lapply(SNPS[cols_snps], function(x) suppressWarnings(as.numeric(x)))

# 3. Filter SNPs (Removed the backslash \$ since this is a pure R file now)
SNPS <- SNPS[SNPS$snp_proportion >= min_snp, ]

# 4. Extract Array Data
Depth <- ASAP.get.depth(ASAP, num_cores = 1)
Proportions <- ASAP.get.proportions(ASAP, num_cores = 1)

# Joining results
array_info <- left_join(Depth, Proportions, by = c("run", "name", "assay_name", "position"))

# 5. Save outputs
save(ASAP, SNPS, array_info, file = paste0(sample_id, "_XML_Data.Rdata"))
write.csv(ASAP, file = paste0(sample_id, "_Summary.csv"), row.names = FALSE)