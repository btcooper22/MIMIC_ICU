require(doParallel)
require(tools)

parallel_setup <- function(n_cores, priority = "idle")
{
  # Set priority
  if(priority == "idle")
  {
    psnice(value = 19)
  } else
  {
    psnice(value = 15)
  }

  registerDoParallel(ifelse(detectCores() <= n_cores,
                            detectCores() - 1,
                            n_cores)
  )
}