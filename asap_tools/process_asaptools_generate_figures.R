#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(plotly)
library(htmlwidgets)
library(zoo) # Required for rolling averages

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: generate_qc_plots.R <rdata> <prefix> <poi_csv>")
}

RDATA_INPUT <- args[1]
PREFIX      <- args[2]
POI_CSV     <- args[3]
SNP_THRESHOLD <- as.numeric(args[4])
SNP_DEPTH     <- as.numeric(args[5])

# RDATA_INPUT <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/f8/c0fa4c75174729c411a3fa2129f452/Combined_ASAP_Data.Rdata"
# PREFIX      <- "Test"
# POI_CSV     <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/tests/Genes_Of_Interest/H37Rv_Genes_Of_Interst.csv"
# SNP_THRESHOLD <- 0.0001
# SNP_DEPTH <- 499
  
load(RDATA_INPUT)

array_info <- final_array

SNPS <- final_snps

# Join metadata
array_info <- left_join(array_info, select(final_asap, run, assay_name, name, amplicon_reads, avg_depth, breadth))

SNPS <- full_join(select(final_array, run, name, assay_name, position, depth), SNPS, by = c("run", "name", "assay_name", "position" = "snp_position")) 
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

# --- Handle Positions of Interest ---
if (!(POI_CSV %in% c("NA", "NULL", "", NA))) {
  genes_poi <- read.csv(POI_CSV)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  genes_expanded <- genes_poi %>%
    rowwise() %>%
    reframe(
      # Keep original metadata
      X = X,
      assay_name = seqnames,
      strand = strand,
      type = type,
      gene = gene,
      # Generate the sequence of positions
      position = get_positions(start, end)
    )
  
  array_info <- genes_expanded %>%
    select(assay_name, position, gene) %>%
    inner_join(array_info, by = c("assay_name", "position")) %>%
    mutate(assay_name = paste(assay_name, "-", gene, sep = ""))
  
  SNPS <- genes_expanded %>%
    select(assay_name, position, gene) %>%
    inner_join(SNPS, by = c("assay_name", "position")) %>%
    mutate(assay_name = paste(assay_name, "-", gene, sep = ""))
  
}


# Dynamic plot data reduction
target_points <- 10000
total_span <- length(unique(array_info$position))
dynamic_k <- max(1, floor(total_span / target_points))


# --- Create Averaged Dataset for Interactive Plots ---
# This calculates the mean of every 10 base pairs to shrink the object size
array_avg <- array_info %>%
  group_by(name, assay_name) %>%
  arrange(position) %>%
  mutate(
    depth_avg = rollmean(depth, k = dynamic_k, fill = NA),
    n_reads_prop = rollmean(100 * (n_reads / depth), k = dynamic_k, fill = NA)
  ) %>%
  filter(!is.na(depth_avg)) %>%
  # Use the dynamic_k to slice the data
  filter(row_number() %% dynamic_k == 0)

# --- Plot 1: Coverage Depth (Interactive) ---
p_cov <- array_avg %>% 
  ggplot(aes(x = position, y = depth_avg, col = name, group = name,
             text = paste0("Sample: ", name,
                           "<br>~Position: ", position,
                           "<br>Mean Depth (10bp): ", round(depth_avg, 1)))) +
    geom_line(alpha = 0.7) + 
    geom_hline(yintercept = SNP_DEPTH, col = "Red", lty = "dashed", alpha = 0.6)+
    facet_wrap(~assay_name, scales = "free") +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = "Reference Depth of Coverage",
      subtitle = paste0("Note: Red dashed line indicates minimum depth (", SNP_DEPTH, "x) for SNP calls."),
      y = paste0("Mean Coverage Depth (", dynamic_k, "bp window, log10 scaled)"),
      x = "Reference Position (BP)"
    )

p_cov

# Save Static (Full Data)
ggsave(paste0(PREFIX, "_coverage_depth.jpg"), width = 12, height = 8, dpi = 300)

# Save Interactive (Averaged Data)
interactive_plot_coverage <- ggplotly(p_cov, tooltip = "text") %>% partial_bundle()
saveWidget(interactive_plot_coverage, paste0(PREFIX, "_coverage_depth.html"), selfcontained = TRUE)

# --- Plot 2: N Read Proportion (Interactive) ---
p_n <- array_avg %>% 
  ggplot(aes(x = position, y = n_reads_prop, col = name, group = name,
             text = paste0("Sample: ", name,
                           "<br>~Position: ", position,
                           "<br>Proporion 'N' Reads (10bp window): ", round(n_reads_prop, 1)))) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~assay_name, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = paste0("Proporion 'N' Reads (", dynamic_k, "bp window, log10 scaled)"),
       x = "Reference Position (BP)",
       title = "Percent Reads with 'N's across Reference",
       subtitle = "Note: 'N's are accumulated from QC, primer masking, and SMOR masking/correction.",)

p_n

# Save Static (Full Data)
ggsave(paste0(PREFIX, "_n_reads_prop.jpg"), width = 12, height = 8, dpi = 300)

# Save Interactive (Averaged Data)
interactive_plot_n_reads <- ggplotly(p_n, tooltip = "text") %>% partial_bundle()
saveWidget(interactive_plot_n_reads, paste0(PREFIX, "_n_reads_prop.html"), selfcontained = TRUE)

# --- Plot 3: SNP Locations (Interactive) ---
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
  mutate(SNP = paste0(snp_reference, position, Call)) %>%
  filter(snp_reference != Call) %>%
  filter(!is.na(snp_proportion))

SNPS <- SNPS %>% 
  filter(depth > SNP_DEPTH)

SNPS$snp_proportion[is.na(SNPS$snp_proportion)] <- 0

SNP_Plot_Data <- SNPS %>% 
  group_by(name, assay_name, position) %>% 
  summarise(Max_SNP_proportion = max(snp_proportion))
  
SNP_Plot_Data_Mean <- SNP_Plot_Data %>% 
  group_by(assay_name, position) %>% 
  summarise(Mean_SNP_proportion = mean(Max_SNP_proportion, na.rm = TRUE), .groups = "drop")

p_SNP <- SNP_Plot_Data %>% 
  ggplot(aes(x = position, y = Max_SNP_proportion, col = name, group = name,
             text = paste0("Sample: ", name,
                           "<br>~Position: ", position,
                           "<br>Max SNP Proportion: ", round(Max_SNP_proportion, 1)))) +
  geom_line(alpha = 0.7) + 
  # Set inherit.aes = FALSE or manually reset col/group to prevent the 'name' error
  geom_line(data = SNP_Plot_Data_Mean, 
            aes(x = position, y = Mean_SNP_proportion), 
            col = "black", 
            lty = "dashed", 
            inherit.aes = FALSE) + 
  facet_wrap(~assay_name, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Max SNP Prevalence (%)",
       x = "Reference Position (BP)",
       title = "SNP prevalence (%) across reference sequences",
       subtitle = "Note: SNPs are expanded and filtered based on user specified SNP thresholds and depths.")

p_SNP

# Save Static (Full Data)
ggsave(paste0(PREFIX, "_SNP_prop.jpg"), width = 12, height = 8, dpi = 300)

# Save Interactive (Averaged Data)
interactive_plot_SNP <- ggplotly(p_SNP, tooltip = "text") %>% partial_bundle()
saveWidget(interactive_plot_SNP, paste0(PREFIX, "_SNP_prop.html"), selfcontained = TRUE)

