#!/usr/bin/env Rscript

library(tidyverse)
library(openxlsx)
library(plotly)
library(htmlwidgets)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: generate_qc_plots.R <rdata> <prefix> <poi_csv>")
}

RDATA_INPUT <- args[1]
PREFIX      <- args[2]
POI_CSV     <- args[3]

load(RDATA_INPUT)
array_info <- final_array

# --- Handle Positions of Interest ---
if (POI_CSV == "NA" || POI_CSV == "NULL" || POI_CSV == "" || is.na(POI_CSV) || is.null(POI_CSV)) {
  positions_of_interest <- unique(array_info$position)
} else {
  genes_poi <- read.csv(POI_CSV)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  positions_of_interest <- unlist(mapply(get_positions, genes_poi$start, genes_poi$end))
  array_info <- array_info %>% filter(position %in% positions_of_interest)
}

# Join metadata
array_info <- left_join(array_info, select(final_asap, run, assay_name, name, amplicon_reads, avg_depth, breadth))

# --- Plot 1: Coverage Depth ---
p <- array_info %>% 
  ggplot(aes(x = position, y = depth, col = name, group = name,
             text = paste0("Sample: ", name,
                           "<br>Position: ", position,
                           "<br>Depth: ", depth,
                           "<br>Amplicon Reads: ", amplicon_reads,
                           "<br>Avg Depth: ", round(avg_depth, 2),
                           "<br>Breadth: ", round(breadth, 2), "%"))) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~assay_name, ncol = 1) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Coverage Depth (X)", x = "Reference Position (BP)") +
  ggtitle("Coverage Depth Across Samples and Assays")

# Save Static JPEG
ggsave(paste0(PREFIX, "_coverage_depth.jpg"), p, width = 12, height = 8, dpi = 300)

# Save Interactive HTML
interactive_plot_coverage <- ggplotly(p, tooltip = "text")
saveWidget(interactive_plot_coverage, paste0(PREFIX, "_coverage_depth.html"), selfcontained = TRUE)

# --- Plot 2: N Reads ---
p1 <- array_info %>% 
  ggplot(aes(x = position, y = n_reads, col = name, group = name,
             text = paste0("Sample: ", name,
                           "<br>Position: ", position,
                           "<br>N_Reads: ", n_reads))) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~assay_name, ncol = 1) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Number of 'N' reads", x = "Reference Position (BP)") +
  ggtitle("'N' Reads Across Samples and Assays")

# Save Static JPEG
ggsave(paste0(PREFIX, "_n_reads.jpg"), p1, width = 12, height = 8, dpi = 300)

# Save Interactive HTML
interactive_plot_n_reads <- ggplotly(p1, tooltip = "text")
saveWidget(interactive_plot_n_reads, paste0(PREFIX, "_n_reads.html"), selfcontained = TRUE)