#' An example function
#'
#' @param optional_argument an argument with a default value
#'
#' @examples
#' # Some example code:
#' df <- example_function()
#' df
#'
#' @export
example_function <- function(optional_argument = 42L){

  # Internal functions are visible within your package code:
  df <- example_internal(optional_argument)
  
  # It is a good idea to use simple S3 classes to allow custom print/summary/autoplot etc methods:
  class(df) <- c("ExampleClass", class(df))

  return(df)
}
