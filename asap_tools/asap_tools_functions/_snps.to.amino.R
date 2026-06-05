snps.to.amino <- function(snp_db, ref_seq, cores = parallelly::availableCores()) {
  library(tidyverse)
  library(genbankr)
  library(Biostrings)
  library(foreach)
  library(doParallel)
  library(parallelly)

  reference    <- suppressWarnings(genbankr::readGenBank(ref_seq))
  Reference_DF <- data.frame(reference@cds)

  # Convert Bioconductor CharacterLists to standard character vectors
  Reference_DF <- as.data.frame(lapply(Reference_DF, function(x) if (is.list(x)) sapply(x, paste, collapse = ";") else x))
  Reference_DF$experiment <- NULL
  Reference_DF$inference  <- NULL
  Reference_DF <- Reference_DF %>%
    mutate(across(where(~inherits(.x, "CharacterList")), ~sapply(.x, paste, collapse = "; ")))

  Reference_DF$sequence <- "No Seq"
  for (i in 1:nrow(Reference_DF)) {
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence),
                                       Reference_DF[i, 2], Reference_DF[i, 3])
  }

  Amino_Acid_List <- snp_db %>%
    separate(SNP, into = c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = F) %>%
    separate(snp_position, into = c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)", remove = F)

  Amino_Acid_List$snp_position <- as.numeric(Amino_Acid_List$snp_position)

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(Amino_Acid_List), .combine = rbind) %dopar% {
    library(dplyr)

    REFERENCE  <- Amino_Acid_List$reference[SNP]
    POSITION   <- Amino_Acid_List$snp_position[SNP]
    MUTATION   <- Amino_Acid_List$snp_mutation[SNP]
    GENOME_SNP <- paste0(REFERENCE, POSITION, MUTATION)

    Out <- data.frame()

    if (MUTATION != "_" & nchar(MUTATION) == 1) {
      for (GENE in 1:nrow(Reference_DF)) {
        if (POSITION >= Reference_DF$start[GENE] & POSITION <= Reference_DF$end[GENE]) {

          SNP_in_gene   <- (POSITION - Reference_DF$start[GENE]) + 1
          Reference_Seq <- Biostrings::DNAString(Reference_DF$sequence[GENE])
          Observed_Seq  <- Biostrings::DNAString(Reference_DF$sequence[GENE])
          Theoretical_Ref <- Reference_Seq[SNP_in_gene]

          Observed_Seq[SNP_in_gene] <- MUTATION
          Observed_SNP <- Observed_Seq[SNP_in_gene]

          if (Reference_DF$strand[GENE] == "-") {
            Reference_Seq   <- Biostrings::reverseComplement(Reference_Seq)
            Observed_Seq    <- Biostrings::reverseComplement(Observed_Seq)
            Theoretical_Ref <- Biostrings::reverseComplement(Theoretical_Ref)
            MUTATION        <- Biostrings::reverseComplement(Observed_SNP)
            SNP_in_gene     <- (Reference_DF$end[GENE] - POSITION) + 1
          }

          AA_Seq      <- Biostrings::translate(Reference_Seq, if.fuzzy.codon = "solve")
          AA_Observed <- Biostrings::translate(Observed_Seq,  if.fuzzy.codon = "solve")

          Mutations <- Biostrings::pairwiseAlignment(AA_Seq, AA_Observed) %>%
            Biostrings::mismatchTable() %>%
            mutate(AA_Change = paste(Reference_DF$gene[GENE], ":", PatternSubstring,
                                     PatternStart, SubjectSubstring, sep = ""))

          Out <- rbind(Out, data.frame(
            SNP                   = as.character(GENOME_SNP),
            snp_position_genome   = POSITION,
            snp_position_gene     = SNP_in_gene,
            snp_mutation          = as.character(MUTATION),
            Theoretical_Reference = as.character(Theoretical_Ref),
            Gene                  = as.character(Reference_DF$gene[GENE]),
            Product               = as.character(Reference_DF$product[GENE]),
            AA                    = as.character(ifelse(
              length(as.character(unique(Mutations$AA_Change))) > 0,
              as.character(unique(Mutations$AA_Change)),
              "Synonymous"
            )),
            SNP_Gene = as.character(paste0(Theoretical_Ref, SNP_in_gene, MUTATION))
          ))
        }
      }
    }
    Out
  }

  stopCluster(cl)

  Out <- full_join(Amino_Acid_List, select(Temp, -snp_mutation))

  Out$snp_mutation <- as.character(unlist(Out$snp_mutation))
  Out$AA[is.na(Out$AA) & nchar(as.character(Out$snp_mutation)) > 1] <- "Insertions Not Supported"
  Out$AA[is.na(Out$AA) & Out$snp_mutation == "_"]                   <- "Deletions Not Supported"
  Out$AA[is.na(Out$AA)]                                             <- "Non-coding SNP"

  Out$SNP_Gene[is.na(Out$SNP_Gene) & nchar(as.character(Out$snp_mutation)) > 1] <- "Insertions Not Supported"
  Out$SNP_Gene[is.na(Out$SNP_Gene) & Out$snp_mutation == "_"]                   <- "Deletions Not Supported"
  Out$SNP_Gene[is.na(Out$SNP_Gene)]                                              <- "Non-coding SNP"

  select(Out, SNP, SNP_Gene, AA, Gene, Product, Theoretical_Reference)
}
