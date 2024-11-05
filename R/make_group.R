#' Generate an object representing a within-group spread model
#'
#' @description
#' This is a wrapper function around the new/initialise methods for R6 and C++
#' implementations of within-group spread models used on the IPDMR course.
#' The method signatures and outputs of the R6 and C++ implementations are
#' (or at least should be!) identical. For more information on available methods
#' see the help file for SEIRclass (the R6 implementation).
#'
#' @param update_type either stochastic or deterministic
#' @param numE the number of sub-compartments desired for the E state
#' @param group_name an optional name for the group
#' @param model_type the compartmental model representation desired (currently only SEIR is supported)
#' @param implementation either C++ or R6
#'
#' @importFrom methods new
#'
#' @examples
#' set.seed(2024-11-05)
#' r6res <- make_group(update_type = "stochastic", implementation = "R6")$run(10, 0.1)
#'
#' set.seed(2024-11-05)
#' cppres <- make_group(update_type = "stochastic", implementation = "C++")$run(10, 0.1)
#'
#' ## Should be identical:
#' stopifnot(unlist(r6res) - unlist(cppres) == 0L)
#'
#' @export
make_group <- function(update_type = c("deterministic","stochastic"), numE = 3L, group_name=NA_character_, model_type = c("SEIR", "SIR", "SI"), implementation = c("C++","R6")){

  update_type <- match.arg(update_type)
  implementation <- toupper(implementation)
  implementation <- match.arg(implementation)
  model_type <- toupper(model_type)
  model_type <- match.arg(model_type)

  stopifnot(model_type=="SEIR")

  if(implementation=="R6"){

    model <- SEIRclass$new(update_type=update_type, numE=numE, group_name=group_name)

  }else if(implementation=="C++"){

    if(update_type=="deterministic"){
      if(numE==3L){
        model <- new(SEIRdet3, 3L, group_name)
      }else{
        model <- new(SEIRdetN, numE, group_name)
      }
    }else{
      if(numE==3L){
        model <- new(SEIRstoc3, 3L, group_name)
      }else{
        model <- new(SEIRstocN, numE, group_name)
      }
    }

  }else{
    stop("Unmatched implementation type")
  }

  return(model)

}
