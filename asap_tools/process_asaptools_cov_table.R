#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: generate_coverage_table.R <combined_rdata> <min_depth> <prefix> <optional_poi_csv>")
}

rdata_input <- args[1]
min_depth   <- as.numeric(args[2])
prefix      <- args[3]
poi_csv     <- args[4] # Will be "NULL" if not provided

poi_csv <- if(length(args) >= 4) args[4] else "NULL"

# rdata_input <- "/scratch/tporter/ASAP_SC2_Validation/ASAP_Illumina_Paired_SE_Test/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "NULL"
# rdata_input <- "/scratch/tporter/ASAP_TB_Validation/ASAP_TB_Subset/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "/scratch/tporter/ASAP_TB_Validation/Updated_TB_Genes.csv"


# 1. Load Data
load(rdata_input) # Loads final_asap, final_snps, final_array

# 2. Handle Positions of Interest (Optional)
if (is.na(poi_csv) || poi_csv == "NULL" || poi_csv == "") {
  message("No Positions of Interest provided. Generating Whole-Reference coverage summary.")

  array_info <- final_array

  Amplicon_Coverage <- array_info %>%
    group_by(name, assay_name) %>%
    summarise(
      total_bp = n(),
      n_cov = sum(depth > min_depth, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, Coverage) %>%
    pivot_wider(names_from = assay_name, values_from = Coverage)


} else {
  message(paste("Loading positions of interest from:", poi_csv))
  genes <- read.csv(poi_csv)

  # Generate positions for each gene in the CSV
  Gene_Positions <- genes %>%
    rowwise() %>%
    do(data.frame(
      position = seq(min(.$start, .$end), max(.$start, .$end)),
      gene = .$gene,
      assay_name = .$seqnames
    )) %>%
    ungroup()

  array_info <- left_join(final_array, Gene_Positions, by = c("position", "assay_name")) %>%
    filter(!is.na(gene)) # Only keep positions that fall within our defined ranges

  Assay_Coverage <- array_info %>%
    group_by(name, assay_name) %>%
    summarise(
      total_bp = n(),
      n_cov = sum(depth > min_depth, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, Coverage) %>%
    pivot_wider(names_from = assay_name, values_from = Coverage)

  POI_Coverage <- array_info %>%
    group_by(name, assay_name, gene) %>%
    summarise(
      total_bp = n(),
      n_cov = sum(depth > min_depth, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, gene, Coverage) %>%
    pivot_wider(names_from = c(assay_name, gene), values_from = Coverage)

  Amplicon_Coverage <- full_join(Assay_Coverage, POI_Coverage)

}

# 5. Styling and Excel Export
wb <- createWorkbook("TGen North")
addWorksheet(wb, "Amplicon_Coverage")

color_breaks_10 <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837")

getStyle_simple_100 <- function(value) {
  # Define the border vector once
  all_borders <- c("top", "bottom", "left", "right")
  if (is.na(value)) {
    return(createStyle(fgFill = "gray75", border = all_borders))
  }
  clamped_value <- max(0, min(100, value))
  color_index <- max(1, min(10, ceiling((clamped_value + 1e-6) / 10)))
  return(createStyle(fgFill = color_breaks_10[color_index], border = all_borders))
}

# Apply styles (Starting from col 2 to skip sample name)
for (row in 1:nrow(Amplicon_Coverage)) {
  for (col in 2:ncol(Amplicon_Coverage)) {
    val <- Amplicon_Coverage[[row, col]]
    addStyle(wb, "Amplicon_Coverage", style = getStyle_simple_100(val), rows = row + 1, cols = col)
  }
}

writeData(wb, "Amplicon_Coverage", Amplicon_Coverage)
saveWorkbook(wb, paste0(prefix, "_Coverage_Report.xlsx"), overwrite = TRUE)
