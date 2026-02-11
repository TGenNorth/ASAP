#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(openxlsx)
library(doParallel)
library(foreach)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: generate_snp_table.R <rdata> <ref>")
}

RDATA_INPUT        <- args[1]
REFERENCE          <- args[2]


# RDATA_INPUT        <- "/scratch/tporter/ASAP_SC2_Validation/ASAP_Illumina_Paired_ASAP_Tools/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# REFERENCE          <- "/tgen_labs/EPIC/tporter/COVID_BAA_202601_Chronic_Bioinformatics/SC2_Reference.gb"


######################
# Load data from ASAP_Import_XML
######################
load(RDATA_INPUT)

SNPS <- final_snps
array_info <- final_array

#################################
## Extract data from files
#################################
# Create cluster
cl <- makeCluster(parallelly::availableCores())
registerDoParallel(cl)

######################
# Extract Unique Amino Acids from SNPs
######################
# Add in distribution for NA values
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0 T=0 C=0 G=0 _=0"

#Find number of SNPS
SNPS <- SNPS %>%
  mutate(space_count = str_count(snp_distribution, " "))

#Extract sub SNPs, and recalculate proportions
SNPS <- SNPS %>%
  relocate(snp_distribution, .after = last_col()) %>%
  separate(snp_distribution, into = paste0("Dist", 1:(1+max(SNPS$space_count))), sep = " ") %>%
  pivot_longer(Dist1:ncol(.), names_to = "Temp", values_to = "Dist") %>%
  select(-Temp) %>%
  filter(!is.na(Dist)) %>%
  separate(Dist, into = c("Call", "n"), sep = "=") %>%
  mutate(snp_proportion = 100*(as.numeric(n)/as.numeric(location_depth))) %>%
  mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
  filter(snp_reference != Call) %>%
  filter(!is.na(snp_proportion))

#Create list of distinct amino acids
SNPS_To_AA <- SNPS %>%
  select(SNP) %>%
  distinct()

# Convert SNP to amino acid
Gene_SNPS <- TGenGenomicTools::genome.snp.to.gene.snp(snp_db = SNPS_To_AA, ref_seq = REFERENCE, cores = parallelly::availableCores())

Gene_SNPS <- Gene_SNPS %>%
  distinct()

Gene_SNPS$Gene_SNP <- Gene_SNPS$SNP_Gene
Gene_SNPS$SNP_Gene <- NULL

Gene_SNPS

Amino_Acids <- TGenGenomicTools::snps.to.amino(snp_db = SNPS_To_AA, ref_seq = REFERENCE, cores = parallelly::availableCores())

Amino_Acids <- select(Amino_Acids, SNP, Product, AA)

save(Amino_Acids, Gene_SNPS, file = paste0("SNP_Amino_Acid_Table.Rdata"))