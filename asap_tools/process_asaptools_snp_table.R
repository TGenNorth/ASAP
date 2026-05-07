#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(openxlsx)
library(doParallel)
library(foreach)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
  stop("Usage: process_asaptools_snp_table.R <rdata> <prefix> <min_snp_perc> <max_snp_count> <min_depth> <remove_names> <poi_csv> <bed_file> <aa_rdata> <ref1> <ref2> ...")
}

RDATA_INPUT        <- args[1]
PREFIX             <- args[2]
MIN_SNP_PERC       <- as.numeric(args[3]) * 100
MAX_SNP_COUNT      <- as.numeric(args[4])
MIN_LOCATION_DEPTH <- as.numeric(args[5])
REMOVE_NAMES       <- if(args[6] == "NONE" || args[6] == "") character(0) else unlist(strsplit(args[6], ","))
POI_CSV            <- args[7]
BED_FILE           <- args[8]
SNP_RDATA          <- args[9]
SNP_XLS            <- args[10]
GB_FILES           <- args[11:length(args)]

# RDATA_INPUT        <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/b4/cee59e06acf4503f68ca757075d169/Combined_ASAP_Data.Rdata"
# PREFIX             <- "RSV_Test"
# MIN_SNP_PERC       <- as.numeric(0.05) * 100
# MAX_SNP_COUNT      <- as.numeric(50)
# MIN_LOCATION_DEPTH <- as.numeric(499)
# REMOVE_NAMES       <- if("NONE" == "NONE" || args[6] == "") character(0) else unlist(strsplit(args[6], ","))
# POI_CSV            <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/tests/Genes_Of_Interest/H37Rv_Genes_Of_Interst.csv"
# BED_FILE           <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/b4/cee59e06acf4503f68ca757075d169/H37Rv_NC0009623_Primer_File_Ampseq.bed"
# SNP_RDATA          <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/b4/cee59e06acf4503f68ca757075d169/SNP_Amino_Acid_Table.Rdata"
# GB_FILES           <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/b4/cee59e06acf4503f68ca757075d169/genbank_input/"

######################
# Load data from ASAP_Import_XML
######################
load(RDATA_INPUT)
SNPS <- final_snps
array_info <- final_array

# --- Handle Flexible Reference List ---
all_gb_paths <- c()

for (path in GB_FILES) {
  if (dir.exists(path)) {
    # If the arg is a directory, get all GenBank files inside
    found_files <- list.files(path, full.names = TRUE, pattern = "\\.(gb|gbk|genbank)$")
    all_gb_paths <- c(all_gb_paths, found_files)
  } else if (file.exists(path)) {
    # If it's a direct file path
    all_gb_paths <- c(all_gb_paths, path)
  }
}

# We use filenames (sans extension) to filter the assay_name later
# e.g., "RSV_A.gb" -> "RSV_A"
valid_refs <- tools::file_path_sans_ext(basename(all_gb_paths))

# Debug message to .command.log so you can verify what was found
if (length(valid_refs) > 0) {
  message(paste("Found", length(valid_refs), "GenBank references:", paste(valid_refs, collapse=", ")))
} else {
  message("No GenBank files found. Proceeding with all SNPs in Rdata.")
}

######################
# Handle Positions of Interest
######################
if (POI_CSV == "NA" || POI_CSV == "NULL" || POI_CSV == "" || is.na(POI_CSV) || is.null(POI_CSV)) {
  message("No Positions of Interest provided. Generating Whole-Genome coverage summary.")
  positions_of_interest <- unique(array_info$position)
} else {
  genes_poi <- read.csv(POI_CSV)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  positions_of_interest <- unlist(mapply(get_positions, genes_poi$start, genes_poi$end))
}

######################
# Load data from bed file (Optional)
######################
if (BED_FILE == "NULL" || !file.exists(BED_FILE) || is.null(BED_FILE)) {
  message("No primer BED file provided. Skipping primer annotation...")
  Primer_Locations <- numeric(0)
} else {
  bed_colnames <- c("reference", "start", "end", "name", "integer", "strand")
  Bed_File <- read_delim(BED_FILE, delim = "\t", col_names = bed_colnames, trim_ws = TRUE)
  get_positions <- function(start, end) seq(start, end)
  Primer_Locations <- unique(unlist(mapply(get_positions, Bed_File$start, Bed_File$end)))
}

######################
# Extract Unique Amino Acids from SNPs
######################
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"
SNPS <- SNPS %>% mutate(space_count = str_count(snp_distribution, " "))
max_spaces <- max(SNPS$space_count, na.rm = TRUE)

SNPS <- SNPS %>%
  relocate(snp_distribution, .after = last_col()) %>%
  separate(snp_distribution, into = paste0("Dist", 1:(1 + max_spaces)), sep = ", ", fill = "right") %>%
  pivot_longer(starts_with("Dist"), names_to = "Temp", values_to = "Dist") %>%
  select(-Temp) %>%
  filter(!is.na(Dist)) %>%
  separate(Dist, into = c("Call", "n"), sep = "=") %>%
  mutate(snp_proportion = 100*(as.numeric(n)/as.numeric(location_depth))) %>%
  mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
  filter(snp_reference != Call) %>%
  filter(!is.na(snp_proportion))

# --- Load Optional AA Data ---
if (SNP_RDATA == "NULL" || !file.exists(SNP_RDATA) || is.null(SNP_RDATA)) {
  message("No SNP Amino Acid Rdata provided. Columns will be empty.")
  Amino_Acids <- data.frame(SNP = character())
  Gene_SNPS   <- data.frame(SNP = character())
} else {
  message("SNP Amino Acid Rdata provided.")
  load(SNP_RDATA) 
  if ("SNP_Gene" %in% colnames(Gene_SNPS)) Gene_SNPS <- rename(Gene_SNPS, Gene_SNP = SNP_Gene)
}

# Join info
AA_Merged <- left_join(Amino_Acids, Gene_SNPS, by = c("assay_name", "SNP"), relationship = "many-to-many")
SNPS <- left_join(SNPS, AA_Merged, by = c("assay_name", "SNP"), relationship = "many-to-many")

######################
# SNP QC & Sample Exclusion
######################
SAMPLE_Exclude <- SNPS %>%
  filter(snp_proportion > MIN_SNP_PERC, location_depth > MIN_LOCATION_DEPTH) %>%
  group_by(name) %>%
  tally() %>%
  filter(n > MAX_SNP_COUNT) %>%
  select(name)

if (length(REMOVE_NAMES) > 0) {
  manual_excl <- SNPS %>%
    filter(grepl(paste(REMOVE_NAMES, collapse = "|"), name)) %>%
    select(name) %>% distinct()
  SAMPLE_Exclude <- rbind(SAMPLE_Exclude, manual_excl) %>% distinct()
}

######################
# Create QCed Table Functions
######################

generate_wide_table <- function(include_only = TRUE) {
  
  # 1. Filter SNPs based on QC and valid references
  sig_filter <- SNPS %>%
    filter(as.numeric(snp_proportion) > MIN_SNP_PERC,
           as.numeric(location_depth) > MIN_LOCATION_DEPTH,
           as.numeric(snp_position) %in% positions_of_interest)
  
  if (length(valid_refs) > 0) {
    sig_filter <- sig_filter %>% 
      filter(grepl(paste(valid_refs, collapse="|"), assay_name))
  }
  
  if (include_only) {
    sig_filter <- sig_filter %>% filter(!name %in% SAMPLE_Exclude$name)
  }
  
  # Get unique significant SNPs for their specific coordinates
  SIGNIFICANT_SNPS_COORDS <- sig_filter %>% 
    select(run, assay_name, name, SNP, snp_position) %>% 
    distinct()
  
  # 2. Prepare Background from array_info
  Background <- array_info %>% distinct()
  
  if (length(valid_refs) > 0) {
    Background <- Background %>% 
      filter(grepl(paste(valid_refs, collapse="|"), assay_name))
  }
  
  # 3. Join by coordinates
  Background <- Background %>%
    left_join(
      SIGNIFICANT_SNPS_COORDS, 
      by = c("run", "assay_name", "name", "position" = "snp_position"),
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(SNP))
  
  # 4. Join back to metadata
  Background <- Background %>%
    left_join(
      SNPS %>% select(run, assay_name, name, SNP, snp_proportion, n, 
                      any_of(c("AA", "Gene_SNP", "Gene", "Product"))),
      by = c("run", "assay_name", "name", "SNP")
    )
  
  # Final cleanup logic
  Background$n[is.na(Background$n)] <- "No SNP Detected"
  Background$snp_proportion[is.na(Background$snp_proportion)] <- 0
  
  # --- Logic for "Unknown" Primer Regions ---
  Background <- Background %>%
    mutate(Primer = case_when(
      # If BED_FILE was NULL or missing, Primer_Locations has length 0
      (BED_FILE == "NULL" || !file.exists(BED_FILE) || is.null(BED_FILE)) ~ "Unknown",
      # Otherwise perform standard T/F check
      position %in% Primer_Locations ~ "TRUE",
      TRUE ~ "FALSE"
    )) %>%
    mutate(snp_prop_final = ifelse(depth > MIN_LOCATION_DEPTH, round(snp_proportion, 2), NA)) %>%
    select(Run = run, 
           Assay = assay_name, 
           Sample = name, 
           any_of("Gene"), 
           `SNP (Genome)` = SNP, 
           `SNP (Gene)` = any_of("Gene_SNP"), 
           `Amino Acid Change` = any_of("AA"), 
           `Primer Region` = Primer, # Will now be TRUE, FALSE, or Unknown
           `SNP Proportion (%)` = snp_prop_final, 
           `SNP Depth` = n, 
           `Location Depth` = depth)
  
  # 5. Pivot to Wide format
  Wide <- Background %>%
    select(-`SNP Depth`, -`Location Depth`) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = `SNP Proportion (%)`)
  
  return(Wide)
}

SNP_TABLE_Included_Samples_Wide <- generate_wide_table(include_only = TRUE)
SNP_Table_Wide_All <- generate_wide_table(include_only = FALSE)

# Save CSVs
write.csv(SNP_TABLE_Included_Samples_Wide, paste0(PREFIX, "_SNP_Table_Included_Samples.csv"))
write.csv(SNP_Table_Wide_All, paste0(PREFIX, "_SNP_Table_All_Samples.csv"))

######################
# Excel Export & Styling
######################

if (SNP_XLS == TRUE){

  cat("Generating excel table, this can take a very long time.")

  getStyle <- function(value) {
    if (is.na(value)) return(createStyle(fgFill = "gray75", border = "TopBottomLeftRight"))
    if (value > 90) return(createStyle(fgFill = "#800026", border = "TopBottomLeftRight"))
    if (value > 50) return(createStyle(fgFill = "#fd8d3c", border = "TopBottomLeftRight"))
    if (value > 0)  return(createStyle(fgFill = "deepskyblue", border = "TopBottomLeftRight"))
    return(createStyle(fgFill = "#9ecae1", border = "TopBottomLeftRight"))
  }

  wb <- createWorkbook()
  sheets <- list("Included_Samples" = SNP_TABLE_Included_Samples_Wide, "All_Samples" = SNP_Table_Wide_All)

  for (sname in names(sheets)) {
    addWorksheet(wb, sname)
    dat <- sheets[[sname]]
    writeData(wb, sname, dat)
    
    # Only apply styling and loops if there are rows present
    if (nrow(dat) > 0) {
      
      # Base Styling
      addStyle(wb, sname, createStyle(border = "TopBottomLeftRight"), 
              rows = 1:(nrow(dat) + 1), cols = 1:ncol(dat), gridExpand = TRUE)
      
      # Header Rotation
      addStyle(wb, sname, createStyle(textRotation = -90), 
              rows = 1, cols = 8:ncol(dat), gridExpand = TRUE)
      
      # Conditional Cell Styling
      for (r in 1:nrow(dat)) {
        for (c in 8:ncol(dat)) {
          addStyle(wb, sname, style = getStyle(dat[r, c]), rows = r + 1, cols = c)
        }
      }
      
      freezePane(wb, sname, firstActiveRow = 2, firstActiveCol = 8)
      
    } else {
      # Optional: Add a message or specific formatting for empty sheets
      message(paste("Sheet", sname, "is empty. Skipping styling loops."))
    }
  }

  saveWorkbook(wb, paste0(PREFIX, "_SNP_Table_Final.xlsx"), overwrite = TRUE)} else {
    cat("Excel file skipped based on input parameters...")
}