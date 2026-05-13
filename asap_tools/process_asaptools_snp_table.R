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

# RDATA_INPUT        <- "/scratch/tporter/ASAP_RSV_Results/work/eb/34919aa89a69685f5f8ec49607d67f/Combined_ASAP_Data.Rdata"
# PREFIX             <- "RSV_Test"
# MIN_SNP_PERC       <- as.numeric(0.9) * 100
# MAX_SNP_COUNT      <- as.numeric(5000000)
# MIN_LOCATION_DEPTH <- as.numeric(499)
# REMOVE_NAMES       <- if("NONE" == "NONE" || args[6] == "") character(0) else unlist(strsplit(args[6], ","))
# POI_CSV            <- NA # "/tgen_labs/EPIC/tporter/ASAP/nextflow/tests/Genes_Of_Interest/H37Rv_Genes_Of_Interst.csv"
# BED_FILE           <- NULL #"/scratch/tporter/ASAP_RSV_Results/work/eb/34919aa89a69685f5f8ec49607d67f/H37Rv_NC0009623_Primer_File_Ampseq.bed"
# SNP_RDATA          <- "/scratch/tporter/ASAP_RSV_Results/work/eb/34919aa89a69685f5f8ec49607d67f/SNP_Amino_Acid_Table.Rdata"
# SNP_XLS            <- FALSE
# GB_FILES           <- "/scratch/tporter/ASAP_RSV_Results/work/eb/34919aa89a69685f5f8ec49607d67f/genbank_input/"

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
if (BED_FILE == "NA" || BED_FILE == "NULL" || is.null(BED_FILE)) { # !file.exists(BED_FILE)
  message("No primer BED file provided. Skipping primer annotation...")
  Primer_Locations <- "No BED file."
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
# Note depending on references there can be 2 records because of overlapping genes.
AA_Merged <- left_join(Amino_Acids, Gene_SNPS, by = c("assay_name", "SNP"), relationship = "many-to-many") %>% 
  distinct()

SNPS <- left_join(SNPS, AA_Merged, by = c("assay_name", "SNP"), relationship = "many-to-many")

######################
# SNP QC & Sample Exclusion
######################
SAMPLE_Exclude <- SNPS %>%
  filter(snp_proportion > MIN_SNP_PERC, location_depth > MIN_LOCATION_DEPTH) %>%
  group_by(assay_name, name) %>%
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
generate_SNP_table <- function(include_only = TRUE) {
  
  # 1. Determine which samples are being processed
  target_samples <- if(include_only) {
    unique(array_info$name[!array_info$name %in% SAMPLE_Exclude$name])
  } else {
    unique(array_info$name)
  }
  
  # 2. Identify "Significant Positions"
  # These are coordinates where at least one sample passed your QC filters
  sig_positions <- SNPS %>%
    filter(as.numeric(snp_proportion) > MIN_SNP_PERC,
           as.numeric(location_depth) > MIN_LOCATION_DEPTH,
           as.numeric(snp_position) %in% positions_of_interest) %>%
    {if (length(valid_refs) > 0) filter(., grepl(paste(valid_refs, collapse="|"), assay_name)) else .} %>%
    select(assay_name, snp_position, SNP, any_of(c("AA", "Gene_SNP", "Gene", "Product"))) %>%
    distinct()
  
  # 3. Join background with SNP metadata.
  Background <- array_info %>%
    filter(name %in% target_samples) %>%
    inner_join(sig_positions, by = c("assay_name", "position" = "snp_position"), relationship = "many-to-many")
  
  # 4. Join with the observed SNP data and SNP calls
  # This tells us if a SPECIFIC sample has that specific SNP.
  Background <- Background %>%
    left_join(
      SNPS %>% select(run, assay_name, name, SNP, snp_proportion, snp_depth),
      by = c("run", "assay_name", "name", "SNP")
    )
  
  Background$snp_proportion[is.na(Background$snp_proportion)] <- 0
  
  # 5. Logical Branching for Coverage vs SNP
  Background <- Background %>%
    mutate(
      snp_prop_final = case_when(
        # CONDITION 1: Depth is too low -> Identify as No Data (NA)
        depth <= MIN_LOCATION_DEPTH ~ paste0("Low Coverage [SNP:", round(snp_proportion, 2), "%, Depth:", depth, "Depth Threshold:", MIN_LOCATION_DEPTH,"]"),
        
        # CONDITION 2: Depth is good and SNP exists or is 0 -> Identify as Variant (%)
        TRUE ~ as.character(round(snp_proportion, 2))
      ),
      Primer = case_when(
        (is.null(BED_FILE) || BED_FILE == "NULL" || BED_FILE == "NA") ~ "No Bed File Provided",
        position %in% Primer_Locations ~ "TRUE",
        TRUE ~ "FALSE"
      )
    )
  
  #6. Pivot to Wide format
  Wide <- Background %>%
    select(Run = run,
           Assay = assay_name,
           Sample = name,
           any_of("Gene"),
           `SNP (Genome)` = SNP,
           `SNP (Gene)` = any_of("Gene_SNP"),
           `Amino Acid Change` = any_of("AA"),
           `Primer Region` = Primer,
           `SNP Proportion (%)` = snp_prop_final) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = `SNP Proportion (%)`)
  
  # 7. Return linelist of SNPS
  SNP_Linelist <- Background %>% 
    filter(depth > MIN_LOCATION_DEPTH) %>% # Reversed logic for clarity: keep if > min
    filter(snp_proportion > MIN_SNP_PERC) %>% 
    filter(SNP %in% sig_positions$SNP)
  
  if (!"Gene" %in% names(Background)) Background$Gene <- "No GB file provided."
  if (!"Gene_SNP" %in% names(Background)) Background$Gene_SNP <- "No GB file provided."
  if (!"AA" %in% names(Background)) Background$AA <- "No GB file provided."
  
  SNP_Linelist <- Background %>% 
    filter(!depth <= MIN_LOCATION_DEPTH) %>% # Filter Low Depth Samples
    filter(snp_proportion > MIN_SNP_PERC) %>% 
    filter(`SNP` %in% sig_positions$SNP) %>% 
    select(run, assay_name, name, `Primer Region` = Primer, `SNP (Genome)` = SNP, Gene, `SNP (Gene)` = `Gene_SNP`, `Amino Acid Change` = AA, `SNP Depth` = snp_depth, `Location Depth` = depth, `SNP Prevalence` = snp_prop_final)

  return(list(Wide, SNP_Linelist))
}

# Generate SNP Tables
SNP_Table_Included <- generate_SNP_table(include_only = TRUE)
SNP_Table_All <- generate_SNP_table(include_only = FALSE)

SNP_Table_Included_Samples_Wide <- SNP_Table_Included[[1]]
SNP_Table_Included_LineList <- SNP_Table_Included[[2]]

SNP_Table_Wide_All <- SNP_Table_All[[1]]
SNP_Table_All_Linelist <- SNP_Table_All[[2]]

# Save CSVs for Wide Data
write.csv(SNP_Table_Included_Samples_Wide, paste0(PREFIX, "_SNP_Table_Included_Samples.csv"))
write.csv(SNP_Table_Wide_All, paste0(PREFIX, "_SNP_Table_All_Samples.csv"))

# Save CSVs for Linelist Data
write.csv(SNP_Table_Included_LineList, paste0(PREFIX, "_SNP_Linelist_Included_Samples.csv"))
write.csv(SNP_Table_All_Linelist, paste0(PREFIX, "_SNP_Linelist_All_Samples.csv"))

######################
# Excel Export & Styling
######################

if (SNP_XLS == TRUE){
  cat("Generating excel table, this can take a very long time.")
  
  getStyle <- function(value) {
    # 1. Handle the "Low Coverage" string check first
    if (grepl("Low Coverage", value)) return(createStyle(fgFill = "gray15", fontColour = "white", border = "TopBottomLeftRight"))
    
    # 2. Convert value to numeric so the '>' comparisons work
    val <- as.numeric(value)
    
    # 3. Use 'val' for comparisons (matching your existing logic)
    if (is.na(val)) return(createStyle(fgFill = "#9ecae1", border = "TopBottomLeftRight")) # Handle non-numeric/NA
    if (val > 99) return(createStyle(fgFill = "#49006a", fontColour = "white", border = "TopBottomLeftRight")) 
    if (val > 90) return(createStyle(fgFill = "#800026", fontColour = "white", border = "TopBottomLeftRight")) 
    if (val > 80) return(createStyle(fgFill = "#bd0026", fontColour = "white", border = "TopBottomLeftRight"))
    if (val > 70) return(createStyle(fgFill = "#e31a1c", fontColour = "white", border = "TopBottomLeftRight"))
    if (val > 50) return(createStyle(fgFill = "#fc4e2a", border = "TopBottomLeftRight")) 
    if (val > 25) return(createStyle(fgFill = "#fd8d3c", border = "TopBottomLeftRight"))
    if (val > 10) return(createStyle(fgFill = "#feb24c", border = "TopBottomLeftRight"))
    if (val > 5)  return(createStyle(fgFill = "#fed976", border = "TopBottomLeftRight")) # Light Yellow-Orange
    if (val > 2)  return(createStyle(fgFill = "#ffeda0", border = "TopBottomLeftRight")) # Pale Yellow
    if (val > 1)  return(createStyle(fgFill = "#ffffcc", border = "TopBottomLeftRight")) # Very Light Cream
    if (val > 0)  return(createStyle(fgFill = "grey75", border = "TopBottomLeftRight")) 
    return(createStyle(fgFill = "#9ecae1", border = "TopBottomLeftRight"))
  }
  
  wb <- createWorkbook()
  sheets <- list("Included_Samples" = SNP_Table_Included_Samples_Wide, "All_Samples" = SNP_Table_Wide_All)
  
  for (sname in names(sheets)) {
    addWorksheet(wb, sname)
    dat <- sheets[[sname]]
    writeData(wb, sname, dat)
    
    if (nrow(dat) > 0) {
      addStyle(wb, sname, createStyle(border = "TopBottomLeftRight"), 
               rows = 1:(nrow(dat) + 1), cols = 1:ncol(dat), gridExpand = TRUE)
      
      addStyle(wb, sname, createStyle(textRotation = -90), 
               rows = 1, cols = 8:ncol(dat), gridExpand = TRUE)
      
      for (r in 1:nrow(dat)) {
        for (c in 8:ncol(dat)) {
          # Use [[ ]] to extract the single value correctly from the data frame
          addStyle(wb, sname, style = getStyle(dat[r, c]), rows = r + 1, cols = c)
        }
      }
      
      freezePane(wb, sname, firstActiveRow = 2, firstActiveCol = 8)
      
    } else {
      message(paste("Sheet", sname, "is empty. Skipping styling loops."))
    }
  }
  
  saveWorkbook(wb, paste0(PREFIX, "_SNP_Table_Final.xlsx"), overwrite = TRUE)
} else {
  cat("Excel file skipped based on input parameters...")
}
