#' Generate an object representing a within-group spread model
#'
#' @description
#' This is a wrapper function around the new/initialise methods for R6 and C++
#' implementations of within-group spread models used on the IPDMR course.
#' For more information on available methods see the documentation for
#' \link{SEIRclass} (the method signatures and outputs of the R6 and C++
#' implementations should be identical).
#'
#' @param update_type either stochastic or deterministic
#' @param numE the number of sub-compartments desired for the E state (0 or more)
#' @param numI the number of sub-compartments desired for the I state (1 or more)
#' @param numR the number of sub-compartments desired for the R state (0 or more)
#' @param group_name an optional name for the group
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
make_group <- function(update_type = c("deterministic","stochastic"), numE = 3L, numI = 1L, numR = 1L, group_name=NA_character_, implementation = c("C++","R6")){

  update_type <- match.arg(update_type)
  implementation <- toupper(implementation)
  implementation <- match.arg(implementation)

  qassert(numE, "X1[0,)")
  qassert(numI, "X1(0,)")
  qassert(numR, "X1[0,)")

  if(implementation=="R6"){

    if(numE==0L) stop("The R6 implementation is limited to numE>=1")
    if(numI!=1L) stop("The R6 implementation is limited to numI==1")
    if(numR!=1L) stop("The R6 implementation is limited to numR==1")

    model <- SEIRclass$new(update_type=update_type, numE=numE, group_name=group_name)

  }else if(implementation=="C++"){

    modname <- str_c("SEIR",
      if(update_type=="deterministic") "det" else "stoc",
      if(numE %in% c(0,1,3)) numE else "N",
      if(numI %in% c(1)) numI else "N",
      if(numR %in% c(0,1)) numR else "N"
    )
    model <- new(get(modname), numE, numI, numR, group_name)

  }else{
    stop("Unmatched implementation type")
  }

  return(model)

}
