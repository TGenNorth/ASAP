read.ASAP.snps.individual <- function(XML) {
  library(xml2)
  library(tidyverse)

  Out      <- data.frame()
  xml_data <- read_xml(XML, options = "HUGE")
  Run_Info <- data.frame(run = "Individual_XML_processing")

  Sample_Node <- xml_data

  Sample_Info <- data.frame(
    name             = xml_attr(Sample_Node, "name"),
    Total_Reads      = xml_attr(Sample_Node, "total_reads"),
    Trimmed_Reads    = xml_attr(Sample_Node, "trimmed_reads"),
    Mapped_Reads     = xml_attr(Sample_Node, "mapped_reads"),
    unassigned_reads = xml_attr(Sample_Node, "unassigned_reads"),
    unmapped_reads   = xml_attr(Sample_Node, "unmapped_reads")
  )

  Assays <- xml_children(Sample_Node)

  for (i in seq_along(Assays)) {
    Assay_Node <- Assays[[i]]

    Assay_Info <- data.frame(
      assay_function = xml_attr(Assay_Node, "function"),
      assay_gene     = xml_attr(Assay_Node, "gene"),
      assay_name     = xml_attr(Assay_Node, "name"),
      assay_type     = xml_attr(Assay_Node, "type")
    )

    Amplicons <- xml_children(Assay_Node)

    for (j in seq_along(Amplicons)) {
      Amplicon_Node <- Amplicons[[j]]
      Amplicon_Info <- data.frame(amplicon_number = j)

      All_Children <- xml_children(Amplicon_Node)
      SNP_Nodes    <- All_Children[xml_name(All_Children) == "snp"]

      if (length(SNP_Nodes) > 0) {
        for (k in seq_along(SNP_Nodes)) {
          SNP_Node <- SNP_Nodes[[k]]

          Dist_Node <- xml_child(SNP_Node, "base_distribution")
          if (!inherits(Dist_Node, "xml_missing")) {
            attrs    <- xml_attrs(Dist_Node)
            SNP_Dist <- paste(names(attrs), attrs, sep = "=", collapse = ", ")
          } else {
            SNP_Dist <- "No reads matching SNP"
          }

          Call_Node <- xml_child(SNP_Node, "snp_call")
          if (!inherits(Call_Node, "xml_missing")) {
            snp_depth      <- xml_attr(Call_Node, "count")
            snp_proportion <- xml_attr(Call_Node, "percent")
            snp_call       <- xml_text(Call_Node)
          } else {
            snp_depth <- "0"; snp_proportion <- "0"; snp_call <- "N/A"
          }

          SNP_Info <- data.frame(
            location_depth   = xml_attr(SNP_Node, "depth"),
            snp_name         = xml_attr(SNP_Node, "name"),
            snp_position     = xml_attr(SNP_Node, "position"),
            snp_reference    = xml_attr(SNP_Node, "reference"),
            snp_depth        = snp_depth,
            snp_proportion   = snp_proportion,
            snp_call         = snp_call,
            snp_distribution = SNP_Dist
          )

          Temp <- cbind(Run_Info, Sample_Info, Assay_Info, Amplicon_Info, SNP_Info)
          Out  <- rbind(Out, Temp)
        }
      }
    }
  }

  print(paste("Processing complete for:", xml_attr(Sample_Node, "name")))

  if (nrow(Out) > 0) {
    row.names(Out) <- 1:nrow(Out)
    num_cols <- c("Total_Reads", "Trimmed_Reads", "Mapped_Reads", "unassigned_reads",
                  "unmapped_reads", "location_depth", "snp_position", "snp_depth", "snp_proportion")
    num_cols <- intersect(num_cols, names(Out))
    Out[num_cols] <- lapply(Out[num_cols], function(x) as.numeric(as.character(x)))
  }

  return(Out)
}
