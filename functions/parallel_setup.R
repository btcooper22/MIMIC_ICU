require(doParallel)
require(tools)

parallel_setup <- function(n_cores)
{
  psnice(value = 19)
  registerDoParallel(ifelse(detectCores() <= n_cores,
                            detectCores() - 1,
                            n_cores)
  )
}