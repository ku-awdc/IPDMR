#' Title
#'
#' @param update_type
#' @param numE
#' @param d_time
#' @param group_name
#' @param implementation
#'
#' @importFrom methods new
#' @export
make_group <- function(update_type = c("deterministic","stochastic"), numE = 3L, d_time=0.1, group_name=NA_character_, model_type = c("SEIR", "SIR", "SI"), implementation = c("R","C++")){

  update_type <- match.arg(update_type)
  implementation <- toupper(implementation)
  implementation <- match.arg(implementation)
  model_type <- toupper(model_type)
  model_type <- match.arg(model_type)

  stopifnot(model_type=="SEIR")

  if(implementation=="R"){

    model <- SEIRclass$new(update_type=update_type, numE=numE, d_time=d_time, group_name=group_name)

  }else if(implementation=="C++"){

    stop("Not implemented yet")
    model <- new(IPDMR:::SEIR,numE,d_time,group_name)

  }else{
    stop("Unmatched implementation type")
  }

  return(model)

}
