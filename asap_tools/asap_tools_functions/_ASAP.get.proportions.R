ASAP.get.proportions <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(foreach)
  library(doParallel)

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  proportions_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run         <- read.ASAP.df$run[[i]]
    name        <- read.ASAP.df$name[[i]]
    assay_name  <- read.ASAP.df$assay_name[[i]]
    proportions <- read.ASAP.df$proportions[[i]]

    proportions <- as.numeric(unlist(strsplit(proportions, ",")))
    position    <- seq_along(proportions)

    data.frame(run, name, assay_name, position, proportions)
  }

  stopCluster(cl)

  proportions_out$position    <- as.numeric(as.character(proportions_out$position))
  proportions_out$proportions <- as.numeric(as.character(proportions_out$proportions))

  return(proportions_out)
}
