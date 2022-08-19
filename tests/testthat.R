library(testthat)
library(recetox.aplcms)

get_num_workers <- function() {
  # CRAN limits the number of cores available to packages to 2
  # source https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions#50571533
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }
  return(num_workers)
}

test_check("recetox.aplcms")
