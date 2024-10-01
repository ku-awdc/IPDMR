## This special R file contains (optional) functions that are run when your package is loaded/attached:

.onAttach <- function(lib, pkg)
{
  ## To print a message uncomment this line:
  # packageStartupMessage("Hello!")
  
  ## To set default settings in the environment we can do e.g.:
  package_env$setting <- 42L
}
