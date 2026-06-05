ASAP.get.nreads <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(foreach)
  library(doParallel)

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  n_reads_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run        <- read.ASAP.df$run[[i]]
    name       <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    n_reads    <- read.ASAP.df$n_reads[[i]]

    n_reads  <- as.numeric(unlist(strsplit(n_reads, ",")))
    position <- seq_along(n_reads)

    data.frame(run, name, assay_name, position, n_reads)
  }

  stopCluster(cl)

  n_reads_out$position <- as.numeric(as.character(n_reads_out$position))

  return(n_reads_out)
}
