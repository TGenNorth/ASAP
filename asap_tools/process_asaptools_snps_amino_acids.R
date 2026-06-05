#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(doParallel)
library(foreach)
library(parallelly)

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_genome.snp.to.gene.snp.R"))
source(file.path(.functions_dir, "_snps.to.amino.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: generate_snp_table.R <rdata> <ref1> <ref2> ...")
}

RDATA_INPUT  <- args[1]
raw_refs <- args[2:length(args)]

# RDATA_INPUT  <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/49/f254b276529b5d69d95d2f6349cb68/Combined_ASAP_Data.Rdata"
# raw_refs <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/49/f254b276529b5d69d95d2f6349cb68/genbank_input/"

GENBANK_FILES <- c()

for (path in raw_refs) {
  if (dir.exists(path)) {
    # If the arg is a directory, get all files inside
    GENBANK_FILES <- c(GENBANK_FILES, list.files(path, full.names = TRUE, pattern = "\\.(gb|gbk|genbank)$"))
  } else if (file.exists(path)) {
    # If it's a direct file path
    GENBANK_FILES <- c(GENBANK_FILES, path)
  }
}

message("Resolved GenBank files:")
print(GENBANK_FILES)

if (length(GENBANK_FILES) == 0) {
  stop("Error: No valid GenBank files found in arguments.")
}

# Load data from ASAP_Import_XML
load(RDATA_INPUT)

SNPS <- final_snps
array_info <- final_array

# --- Initial SNPS Cleaning (as per your original code) ---
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

#Find number of SNPS
SNPS <- SNPS %>%
  mutate(space_count = str_count(snp_distribution, " "))

#Extract sub SNPs, and recalculate proportions
SNPS <- SNPS %>%
  relocate(snp_distribution, .after = last_col()) %>%
  separate(snp_distribution, into = paste0("Dist", 1:(1+max(SNPS$space_count))), sep = ", ") %>%
  pivot_longer(Dist1:ncol(.), names_to = "Temp", values_to = "Dist") %>%
  select(-Temp) %>%
  filter(!is.na(Dist)) %>%
  separate(Dist, into = c("Call", "n"), sep = "=") %>%
  mutate(snp_proportion = 100*(as.numeric(n)/as.numeric(location_depth))) %>%
  mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
  filter(snp_reference != Call) %>%
  filter(!is.na(snp_proportion))

# --- Loop Through All GenBank Files ---
all_amino_acids <- list()
all_gene_snps <- list()

for (REFERENCE in GENBANK_FILES) {
  # Read the specific reference
  reference_obj <- suppressWarnings(genbankr::readGenBank(REFERENCE))

  acc_id <- reference_obj@accession

  file_base <- tools::file_path_sans_ext(basename(REFERENCE))

  message(paste0("Processing: ", REFERENCE, " (ID: ", acc_id, " | FileBase: ", file_base, ")"))

  # Filter SNPs belonging to this specific accession/assay
  SNPS_To_AA <- SNPS %>%
    filter(grepl(acc_id, assay_name) | grepl(file_base, assay_name)) %>%
    select(SNP, assay_name) %>%
    distinct()

  if (nrow(SNPS_To_AA) == 0) {
    message(paste("No SNPs found for accession:", acc_id))
    message(paste("If this is unexpected check SNP assay name."))
    message(paste0("Filtering was conducted with: ", acc_id, " | FileBase: ", file_base, ")"))
    next
  }

  message(paste0("Processing Gene SNPS: ", REFERENCE, " (ID: ", acc_id, " | FileBase: ", file_base, ")"))

  # Convert SNP to Gene SNP using the specific reference
  gene_snps_sub <- suppressWarnings(genome.snp.to.gene.snp(
    snp_db = SNPS_To_AA,
    ref_seq = REFERENCE,
    cores = parallelly::availableCores()
  ) %>%
    mutate(assay_name=file_base))

  message(paste0("Processing Amino Acids: ", REFERENCE, " (ID: ", acc_id, " | FileBase: ", file_base, ")"))

  # Convert SNP to amino acid using the specific reference
  amino_acids_sub <- suppressWarnings(snps.to.amino(
    snp_db = SNPS_To_AA,
    ref_seq = REFERENCE,
    cores = parallelly::availableCores()
  ) %>%
    mutate(assay_name=file_base))

  all_gene_snps[[acc_id]] <- gene_snps_sub
  all_amino_acids[[acc_id]] <- amino_acids_sub

  message(paste0("Success processing: ", REFERENCE, " (ID: ", acc_id, " | FileBase: ", file_base, ")"))
}

# Combine results
Gene_SNPS   <- bind_rows(all_gene_snps) %>% distinct()
Amino_Acids <- bind_rows(all_amino_acids) %>% select(assay_name, SNP, Product, AA)

save(Amino_Acids, Gene_SNPS, file = "SNP_Amino_Acid_Table.Rdata")
