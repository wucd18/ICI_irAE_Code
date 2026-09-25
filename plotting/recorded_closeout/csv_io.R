# CSV data carry a UTF-8 BOM for spreadsheet compatibility on Windows.
# Read with an explicit encoding so the first field never acquires a BOM prefix.
read.csv <- function(file, ...) {
  args <- list(...)
  normal_names <- if(is.null(args$check.names)) TRUE else args$check.names
  args$check.names <- FALSE
  args$file <- file
  args$encoding <- 'UTF-8'
  x <- do.call(utils::read.csv, args)
  names(x) <- sub('^\ufeff', '', names(x))
  if(normal_names) names(x) <- make.names(names(x),unique=TRUE)
  x
}
