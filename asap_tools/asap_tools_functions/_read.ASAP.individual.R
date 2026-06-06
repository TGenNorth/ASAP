read.ASAP.individual <- function(XML) {
  library(xml2)
  library(tidyverse)

  Out <- data.frame()
  xml_data <- xml2::read_xml(XML, options = "HUGE")

  Run_Info <- data.frame(run = "Individual_XML_processing")

  Sample_Node <- xml_data

  Sample_Info <- data.frame(
    name             = xml_attr(Sample_Node, "name"),
    total_reads      = xml_attr(Sample_Node, "total_reads"),
    trimmed_reads    = xml_attr(Sample_Node, "trimmed_reads"),
    mapped_reads     = xml_attr(Sample_Node, "mapped_reads"),
    unassigned_reads = xml_attr(Sample_Node, "unassigned_reads"),
    unmapped_reads   = xml_attr(Sample_Node, "unmapped_reads")
  )

  Assays <- xml_children(Sample_Node)
  if (length(Assays) == 0) return(NULL)

  for (i in 1:length(Assays)) {
    Assay_Node <- Assays[[i]]

    Assay_Info <- data.frame(
      assay_function = xml_attr(Assay_Node, "function"),
      assay_gene     = xml_attr(Assay_Node, "gene"),
      assay_name     = xml_attr(Assay_Node, "name"),
      assay_type     = xml_attr(Assay_Node, "type")
    )

    Amplicons <- xml_children(Assay_Node)
    if (length(Amplicons) == 0) next

    for (AMPLICON in 1:length(Amplicons)) {
      Amplicon_Node <- Amplicons[[AMPLICON]]

      get_node_text <- function(parent, child_name, default) {
        node <- xml_child(parent, child_name)
        if (inherits(node, "xml_missing")) return(default)
        return(as.character(xml_contents(node)))
      }

      Amplicon_Info <- data.frame(
        amplicon_number      = AMPLICON,
        amplicon_reads       = xml_attr(Amplicon_Node, "reads"),
        aligned_reads        = xml_attr(Amplicon_Node, "aligned_reads"),
        primer_reads         = xml_attr(Amplicon_Node, "primer_reads"),
        no_primer_reads      = xml_attr(Amplicon_Node, "no_primer_reads"),
        identity_input       = xml_attr(Amplicon_Node, "identity_input"),
        identity_discarded   = xml_attr(Amplicon_Node, "identity_discarded"),
        smor_input           = xml_attr(Amplicon_Node, "smor_input"),
        smor_pairs_dropped   = xml_attr(Amplicon_Node, "smor_pairs_dropped"),
        smor_consensus_reads = xml_attr(Amplicon_Node, "smor_consensus_reads"),
        amplicon_variant = ifelse(!is.na(xml_attr(Amplicon_Node, "variant")), xml_attr(Amplicon_Node, "variant"), "No variant"),
        breadth          = get_node_text(Amplicon_Node, "breadth", "No Breadth"),
        avg_depth        = get_node_text(Amplicon_Node, "average_depth", "No Average Depth"),
        consensus_seq    = get_node_text(Amplicon_Node, "consensus_sequence", "No Consensus Sequence"),
        depths           = get_node_text(Amplicon_Node, "depths", "No Depth"),
        proportions      = get_node_text(Amplicon_Node, "proportions", "No Proportions"),
        quality_discards = get_node_text(Amplicon_Node, "quality_discards", "No QC Analysis"),
        n_reads          = get_node_text(Amplicon_Node, "n_reads", "No QC Analysis")
      )

      Temp <- cbind(Run_Info, Sample_Info, Assay_Info, Amplicon_Info)
      Out  <- rbind(Out, Temp)
    }

    print(paste("Processing complete for assay:", xml_attr(Assay_Node, "name")))
  }

  num_cols <- c("total_reads", "trimmed_reads", "mapped_reads", "unassigned_reads",
                "unmapped_reads", "amplicon_reads", "aligned_reads",
                "primer_reads", "no_primer_reads",
                "identity_input", "identity_discarded",
                "smor_input", "smor_pairs_dropped", "smor_consensus_reads",
                "breadth", "avg_depth")
  num_cols <- intersect(num_cols, names(Out))
  Out[num_cols] <- lapply(Out[num_cols], function(x) as.numeric(as.character(x)))

  return(Out)
}
