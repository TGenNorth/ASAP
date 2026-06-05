ASAP.get.quality.discards <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(foreach)
  library(doParallel)

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  quality_discards_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run             <- read.ASAP.df$run[[i]]
    name            <- read.ASAP.df$name[[i]]
    assay_name      <- read.ASAP.df$assay_name[[i]]
    quality_discards <- read.ASAP.df$quality_discards[[i]]

    quality_discards <- as.numeric(unlist(strsplit(quality_discards, ",")))
    position         <- seq_along(quality_discards)

    data.frame(run, name, assay_name, position, quality_discards)
  }

  stopCluster(cl)

  quality_discards_out$position         <- as.numeric(as.character(quality_discards_out$position))
  quality_discards_out$quality_discards <- as.numeric(as.character(quality_discards_out$quality_discards))

  return(quality_discards_out)
}
