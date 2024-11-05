#' Combine compartment sizes and rates to obtain compartment change values
#'
#' @description
#' This is utility function taking a vector of compartments and rates, and
#' returning a matrix of calculated (for deterministic) or sampled (for
#' stochastic) values for the change in compartment size from each input
#' compartment (rows) to each of the output compartments represented by the
#' input rates (columns). For most cases there will be only 1 compartment,
#' so a matrix with a single row is returned. Multiple rows are only relevant
#' where there are multiple sub-compartments, i.e. movement from the E state.
#'
#' @param compartments a vector of compartment sizes on which to apply all rates
#' @param rates a vector of rates representing movement from each of the compartments
#' @param d_time the size of the time step
#' @param update_type deterministic or stochastic
#'
#' @export
apply_rates <- function(compartments, rates, d_time, update_type=c("deterministic","stochastic")){

  qassert(compartments, "N+[0,)")
  qassert(rates, "N+[0,)")
  qassert(d_time, "N1(0,)")
  update_type <- match.arg(update_type)

  ## Generate proportions and remove overlap:
  leave <- 1-exp(-sum(rates)*d_time)
  if(leave==0){
    prop <- rep(leave, length(rates))
  }else{
    prop <- leave*rates/sum(rates)
  }

  ## Convert to numbers
  if(update_type=="stochastic"){
    if(!all(compartments%%1==0)) stop("One or more input compartments was not discrete")

    prob <- c(1-sum(prop), prop)
    vapply(compartments, \(x){
      rmultinom(1, x, prob)[-1L,,drop=TRUE]
    }, numeric(length(prop))) |>
      t() ->
      rv
    stopifnot(apply(rv,1,sum)<=compartments)

  }else if(update_type=="deterministic"){

    nb <- length(compartments)
    size <- matrix(rep(compartments, times=length(prop)), nrow=nb)
    prob <- matrix(rep(prop, each=nb), nrow=nb)
    rv <- size*prob

    stopifnot(apply(rv,1,sum) <= compartments)
  }else{
    stop("Unrecognised update_type")
  }

  stopifnot(dim(rv)==c(length(compartments),length(rates)))

  if(!is.null(names(rates))){
    colnames(rv) <- names(rates)
  }

  return(rv)
}
