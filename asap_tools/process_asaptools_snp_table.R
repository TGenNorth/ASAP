#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(openxlsx)
library(doParallel)
library(foreach)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop("Usage: generate_snp_table.R <rdata> <prefix> <min_snp_perc> <max_snp_count> <min_depth> <remove_names> <poi_csv> <ref> <primer_bed>")
}

RDATA_INPUT        <- args[1]
PREFIX             <- args[2]
MIN_SNP_PERC       <- as.numeric(args[3])
MAX_SNP_COUNT      <- as.numeric(args[4])
MIN_LOCATION_DEPTH <- as.numeric(args[5])
REMOVE_NAMES       <- if(args[6] == "NONE" || args[6] == "") character(0) else unlist(strsplit(args[6], ","))
POI_CSV            <- args[7]
REFERENCE          <- args[8]
BED_FILE           <- args[9]


# RDATA_INPUT        <- "/scratch/tporter/ASAP_SC2_Validation/ASAP_Illumina_Paired_ASAP_Tools/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# PREFIX             <- "ASAP_Illumina_Paired_ASAP_Tools"
# MIN_SNP_PERC       <- 0.01
# MAX_SNP_COUNT      <- 5000
# MIN_LOCATION_DEPTH <- 99
# REMOVE_NAMES       <- character(0)
# POI_CSV            <- NA
# REFERENCE          <- "/tgen_labs/EPIC/tporter/COVID_BAA_202601_Chronic_Bioinformatics/SC2_Reference.gb"
# BED_FILE           <- "/tgen_labs/EPIC/tporter/COVID_BAA_202601_Chronic_Bioinformatics/bedfiles/Illumina/Concatenated_Bedfile.bed"

######################
# Load data from ASAP_Import_XML
######################
load(RDATA_INPUT)

SNPS <- final_snps
array_info <- final_array

if (is.na(POI_CSV) || POI_CSV == "NULL" || POI_CSV == "") {
  message("No Positions of Interest provided. Generating Whole-Genome coverage summary.")
  
  # Create a Gene_Positions table for the entire range of the reference
  # We look at final_array to see what positions were actually captured
  genome_min <- min(final_array$position, na.rm = TRUE)
  genome_max <- max(final_array$position, na.rm = TRUE)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  positions_of_interest <- unlist(mapply(get_positions, genome_min, genome_max))
  
} else {
  genes_poi <- read.csv(POI_CSV)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  positions_of_interest <- unlist(mapply(get_positions, genes_poi$start, genes_poi$end))
}

######################
# Load data from bed file and get primer locations
######################
# Function to generate positions for each gene
get_positions <- function(start, end) {
  return(seq(start, end))
}
# Define the column names as requested
bed_colnames <- c("reference", "start", "end", "name", "integer", "strand")

# Read the file
Bed_File <- read_delim(
  BED_FILE, 
  delim = "\t", 
  col_names = bed_colnames,
  trim_ws = TRUE
)

Primer_Locations <-  unique(unlist(mapply(get_positions, Bed_File$start, Bed_File$end)))

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

Amino_Acids <- left_join(Amino_Acids, Gene_SNPS,  relationship = "many-to-many")

SNPS <- left_join(SNPS, Amino_Acids,  relationship = "many-to-many")

######################
# SNP QC
######################
SNPS %>%
  filter(snp_proportion > MIN_SNP_PERC) %>%
  filter(location_depth > MIN_LOCATION_DEPTH) %>%
  group_by(name) %>%
  tally() %>%
  arrange(desc(n)) %>%
  ggplot(aes(x = n))+
  geom_histogram(bins = 100) +
  theme_bw()+
  xlab("Number of SNPs per sample")

SAMPLE_Exclude <- SNPS %>%
  filter(snp_proportion > MIN_SNP_PERC) %>%
  filter(location_depth > MIN_LOCATION_DEPTH) %>%
  group_by(name) %>%
  tally() %>%
  filter(n > MAX_SNP_COUNT) %>%
  select(name)

SAMPLE_Exclude

SAMPLE_Exclude <- rbind(SAMPLE_Exclude,
                        SNPS %>%
                          {if (length(REMOVE_NAMES) > 0) 
                            filter(., grepl(paste(REMOVE_NAMES, collapse = "|"), name)) 
                            else filter(., FALSE)} %>%  # Returns the structure but 0 rows
                          select(name) %>%
                          distinct())

######################
# Define Color Functions
######################

# Function to get the appropriate style based on the value
getStyle <- function(value) {
  if (is.na(value)) {
    return(createStyle(fgFill = "gray75", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 9*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#800026", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 8*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#bd0026", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 7*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#e31a1c", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 6*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#fc4e2a", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 5*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#fd8d3c", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 4*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#feb24c", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 3*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#fed976", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 2*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#ffeda0", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 1*MIN_SNP_PERC) {
    return(createStyle(fgFill = "#ffffcc", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else if (value > 0) {
    return(createStyle(fgFill = "deepskyblue", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  } else {
    return(createStyle(fgFill = "#9ecae1", border = c("top", "bottom", "left", "right"), borderColour = "black"))
  }
}

######################
# Create QCed Table
######################

SIGNIFICANT_SNPS <- SNPS %>%
  filter(!name %in% SAMPLE_Exclude$name) %>%
  filter(as.numeric(snp_proportion) > MIN_SNP_PERC) %>%
  filter(as.numeric(location_depth) > MIN_LOCATION_DEPTH) %>%
  filter(as.numeric(snp_position) %in% positions_of_interest) %>%
  select(SNP) %>%
  distinct()

Background <- full_join(array_info %>% 
                          distinct(),
                        Amino_Acids %>%
                          separate(SNP, into =c("reference", "position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
                          separate(position, into =c("position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F) %>%
                          mutate(position = as.numeric(position)))

Background <- left_join(Background,
                        SNPS %>%
                          select(run, name, SNP, snp_proportion, n, location_depth) %>%
                          filter(SNP %in% SIGNIFICANT_SNPS$SNP),
                        by = c("run", "name", "SNP"))

Background <- Background %>%
  filter(!is.na(SNP))

Background$n[is.na(Background$n)] <-"No SNP Detected"

Background$snp_proportion[is.na(Background$snp_proportion)] <- 0

Background <- Background %>%
  mutate(Primer = ifelse(position %in% Primer_Locations, TRUE, FALSE)) %>%
  select(run, name, Gene, SNP, Gene_SNP, AA, Primer, snp_proportion, n, depth) %>%
  mutate(snp_proportion = round(snp_proportion,2))

names(Background) <- c("Run", "Sample", "Gene", "SNP (Genome)", "SNP (Gene)", "Amino Acid Change", "Primer Region", "SNP Proportion (%)", "SNP Depth", "Location Depth")

Background <- Background %>%
  mutate(`SNP Proportion (%)` = ifelse(`Location Depth` > MIN_LOCATION_DEPTH, `SNP Proportion (%)`, NA))

SNP_TABLE_Included_Samples_Wide <- Background %>%
  filter(!Sample %in% SAMPLE_Exclude$name)%>%
  filter(`SNP (Genome)` %in% SIGNIFICANT_SNPS$SNP) %>%
  select(-`SNP Depth`,-`Location Depth`) %>% 
  distinct() %>% 
  pivot_wider(names_from = Sample, values_from = `SNP Proportion (%)`)

library(openxlsx)

wb <- createWorkbook("ASAP Output")
addWorksheet(wb, "SNP_TABLE_Included_Samples")
writeData(wb, "SNP_TABLE_Included_Samples", SNP_TABLE_Included_Samples_Wide)

# Apply the border style to all cells in the data range
addStyle(wb, sheet = "SNP_TABLE_Included_Samples", style = createStyle(border = c("top", "bottom", "left", "right"), borderColour = "black"), rows = 1:(nrow(SNP_TABLE_Included_Samples_Wide) + 1), cols = 1:ncol(SNP_TABLE_Included_Samples_Wide), gridExpand = TRUE)

# Apply the rotation style to the first row
addStyle(wb, sheet = "SNP_TABLE_Included_Samples", style = createStyle(textRotation = -90, border = c("top", "bottom", "left", "right"), borderColour = "black"), rows = 1, cols = 6:ncol(SNP_TABLE_Included_Samples_Wide), gridExpand = TRUE)

for (row in 1:nrow(SNP_TABLE_Included_Samples_Wide)) {
  for (col in 7:ncol(SNP_TABLE_Included_Samples_Wide)) {
    value <- SNP_TABLE_Included_Samples_Wide[row, col]
    style <- getStyle(value)
    addStyle(wb, sheet = "SNP_TABLE_Included_Samples", style = style, rows = row + 1, cols = col, gridExpand = TRUE)
  }
}

# Freeze panes: freeze between row 1 and 2 and columns 5 and 6
freezePane(wb, sheet = "SNP_TABLE_Included_Samples", firstActiveRow = 2, firstActiveCol = 7)

######################
# Create All Sample Table
######################

SIGNIFICANT_SNPS <- SNPS %>%
  #filter(name %in% SAMPLE_Exclude$name) %>%
  filter(as.numeric(snp_proportion) > MIN_SNP_PERC) %>%
  filter(as.numeric(location_depth) > MIN_LOCATION_DEPTH) %>%
  filter(as.numeric(snp_position) %in% positions_of_interest) %>%
  select(SNP) %>%
  distinct()

Background <- full_join(array_info,
                        Amino_Acids %>%
                          separate(SNP, into =c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
                          separate(snp_position, into =c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F) %>%
                          mutate(snp_position = as.numeric(snp_position)),
                        by = join_by("position" == "snp_position"))

Background

SNPS

Background <- left_join(Background,
                        SNPS %>%
                          select(run, name, SNP, snp_proportion, n, location_depth) %>%
                          filter(SNP %in% SIGNIFICANT_SNPS$SNP),
                        by = c("run", "name", "SNP"))

Background <- Background %>%
  filter(!is.na(SNP))

Background$n[is.na(Background$n)] <-"No SNP Detected"

Background$snp_proportion[is.na(Background$snp_proportion)] <- 0

Background <- Background %>%
  mutate(Primer = ifelse(position %in% Primer_Locations, TRUE, FALSE)) %>%
  select(run, name, Gene, SNP, Gene_SNP, AA, Primer, snp_proportion, n, depth) %>%
  mutate(snp_proportion = round(snp_proportion,2))

names(Background) <- c("Run", "Sample", "Gene", "SNP (Genome)", "SNP (Gene)", "Amino Acid Change", "Primer Region", "SNP Proportion (%)", "SNP Depth", "Location Depth")

unique(Background$Sample)

Background <- Background %>%
  mutate(`SNP Proportion (%)` = ifelse(`Location Depth` > MIN_LOCATION_DEPTH, `SNP Proportion (%)`, NA))

SNP_Table_Wide_All <- Background %>%
  filter(`SNP (Genome)` %in% SIGNIFICANT_SNPS$SNP) %>%
  select(-`SNP Depth`,-`Location Depth`) %>%
  distinct() %>%
  pivot_wider(names_from = Sample, values_from = `SNP Proportion (%)`)

addWorksheet(wb, "SNP_Table_All_Samples")
writeData(wb, "SNP_Table_All_Samples", SNP_Table_Wide_All)

# Apply the border style to all cells in the data range
addStyle(wb, sheet = "SNP_Table_All_Samples", style = createStyle(border = c("top", "bottom", "left", "right"), borderColour = "black"), rows = 1:(nrow(SNP_Table_Wide_All) + 1), cols = 1:ncol(SNP_Table_Wide_All), gridExpand = TRUE)

# Apply the rotation style to the first row
addStyle(wb, sheet = "SNP_Table_All_Samples", style = createStyle(textRotation = -90, border = c("top", "bottom", "left", "right"), borderColour = "black"), rows = 1, cols = 6:ncol(SNP_Table_Wide_All), gridExpand = TRUE)

for (row in 1:nrow(SNP_Table_Wide_All)) {
  for (col in 7:ncol(SNP_Table_Wide_All)) {
    value <- SNP_Table_Wide_All[row, col]
    style <- getStyle(value)
    addStyle(wb, sheet = "SNP_Table_All_Samples", style = style, rows = row + 1, cols = col, gridExpand = TRUE)
  }
}

# Freeze panes: freeze between row 1 and 2 and columns 5 and 6
freezePane(wb, sheet = "SNP_Table_All_Samples", firstActiveRow = 2, firstActiveCol = 7)

saveWorkbook(wb,  paste0(PREFIX, "_SNP_Table_Export_MINSNP_", MIN_SNP_PERC,"_MAXSNPCOUNT_", MAX_SNP_COUNT,"_MINDEPTH_",MIN_LOCATION_DEPTH,".xlsx"), overwrite  = TRUE)

