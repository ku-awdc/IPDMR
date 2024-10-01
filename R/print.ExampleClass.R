## This is an example of a simple print method for an S3 class:

#' @export
print.ExampleClass <- function(x, ...){
  
  # Defensive programming:
  stopifnot("Incomplete ExampleClass object" = "Value" %in% names(x))
  
  # Specific things for this class:
  cat("Hello from an ExampleClass with value ", x[["Value"]], "!\n\n", sep="")
  
  # Pass through to the print method for data frames:
  print.data.frame(x)
  
  # It is recommended to return the object invisibly:
  invisible(x)
  
}
