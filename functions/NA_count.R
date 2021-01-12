require(foreach)

NA_count <- function(.df, start = 2)
{
  foreach(i = 2:ncol(.df), .combine = "rbind") %do%
    {
      data.frame(name = names(.df)[i], 
                 NA_count = sum(is.na(.df[,i])))
    }
}

