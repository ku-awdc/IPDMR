## This function is not exported i.e. it is internal to your package

## We use functions from other packages, so we must import them
## AND also include them in the DESCRIPTION under Imports:
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom rlang .data

example_internal <- function(argument){

  ## Check that the argument is numeric and not missing:
  stopifnot(is.numeric(argument), !is.na(argument))

  ## It is good practice to use the .data pronoun within dplyr functions:
  df <- tibble(Variable = "Test", Value = argument) |>
    mutate(Value2 = .data$Value^2)

  return(df)

}


## If we need the package to remember things within this R session, we can set them in an environment:
package_env <- new.env()
