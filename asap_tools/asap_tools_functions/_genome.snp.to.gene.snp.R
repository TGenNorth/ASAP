genome.snp.to.gene.snp <- function(snp_db, ref_seq, cores = parallelly::availableCores()) {
  library(tidyverse)
  library(genbankr)
  library(Biostrings)
  library(foreach)
  library(doParallel)
  library(parallelly)

  reference    <- suppressWarnings(genbankr::readGenBank(ref_seq))
  Reference_DF <- left_join(
    data.frame(reference@genes),
    data.frame(reference@cds) %>% select(locus_tag, product, translation)
  )

  Reference_DF$sequence <- "No Seq"
  for (i in 1:nrow(Reference_DF)) {
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence),
                                       Reference_DF[i, 2], Reference_DF[i, 3])
  }

  Reference_DF <- Reference_DF %>%
    mutate(gene = ifelse(is.na(gene), locus_tag, gene))

  SNP_List <- snp_db %>%
    separate(SNP, into = c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
    separate(snp_position, into = c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F)

  SNP_List$snp_position <- as.numeric(SNP_List$snp_position)

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(SNP_List), .combine = rbind) %dopar% {
    library(dplyr)

    REFERENCE  <- SNP_List$reference[SNP]
    POSITION   <- SNP_List$snp_position[SNP]
    MUTATION   <- SNP_List$snp_mutation[SNP]
    GENOME_SNP <- paste0(REFERENCE, POSITION, MUTATION)

    Out <- data.frame()

    for (GENE in 1:nrow(Reference_DF)) {
      if (POSITION >= Reference_DF$start[GENE] & POSITION <= Reference_DF$end[GENE]) {

        SNP_in_gene   <- (POSITION - Reference_DF$start[GENE]) + 1
        Reference_Seq <- Biostrings::DNAString(Reference_DF$sequence[GENE])

        if (MUTATION != "_") {
          Observed_Seq <- Biostrings::DNAString(MUTATION)
        } else {
          Observed_Seq <- "_"
        }

        Theoretical_Ref <- Reference_Seq[SNP_in_gene]

        if (Reference_DF$strand[GENE] == "-") {
          SNP_in_gene     <- (Reference_DF$end[GENE] - POSITION) + 1
          Reference_Seq   <- Biostrings::reverseComplement(Reference_Seq)
          Theoretical_Ref <- Biostrings::reverseComplement(Theoretical_Ref)

          if (MUTATION != "_") {
            Observed_Seq <- Biostrings::reverseComplement(Observed_Seq)
          } else {
            Observed_Seq <- "_"
          }

          MUTATION <- Observed_Seq
        }

        Out <- rbind(Out, data.frame(
          SNP                 = GENOME_SNP,
          snp_position_genome = POSITION,
          snp_position_gene   = SNP_in_gene,
          snp_mutation        = MUTATION,
          Theoretical_Reference = Theoretical_Ref,
          Gene                = Reference_DF$gene[GENE],
          SNP_Gene            = paste0(Theoretical_Ref, SNP_in_gene, MUTATION)
        ))
      }
    }
    Out
  }

  stopCluster(cl)

  Out <- full_join(select(SNP_List, SNP), Temp)
  Out$Gene[is.na(Out$Gene)]         <- "Non-gene region"
  Out$SNP_Gene[is.na(Out$SNP_Gene)] <- "Non-gene region"

  select(Out, SNP, Gene, SNP_Gene)
}
