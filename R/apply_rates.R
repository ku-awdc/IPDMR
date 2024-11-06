#' Combine compartment sizes and rates to obtain compartment change values
#'
#' @description
#' This is utility function taking a vector of compartments and rates, and
#' returning a matrix of calculated (for deterministic) or sampled (for
#' stochastic) values for the change in compartment size from each input
#' compartment (row) to each of the output compartments represented by the
#' input rates (columns). When only a single rate is input then the function
#' just returns the usual compartment * 1-exp(-rate*d_time) (for deterministic) or
#' rbinom(length(compartment), compartment, 1-exp(-rate*d_time)) (for stochastic).
#' Where there are multiple rates, then the function corrects for the overlap
#' in the competing rates and (for stochastic) uses a multinomial distribution
#' instead of a binomial distribution to ensure that the sum of each row of
#' output is less than or equal to the corresponding input compartment size.
#'
#' Multiple rows are only relevant where there are multiple sub-compartments,
#' i.e. movement from the E state. In all other cases, a matrix with a single
#' row is returned. Likewise, multiple columns are only relevant where there are
#' multiple rates - in all other cases, a matrix with a single column is returned.
#' If you know that you have a single compartment and/or a single rate, then
#' you can either explicitly select the first row and/or column or use as.numeric()
#' on the output to remove the dimensions attribute (see the bottom of the
#' examples).
#'
#' @param compartments a vector of compartment sizes on which to apply all rates
#' @param rates a vector of rates representing movement from each of the compartments
#' @param d_time the size of the time step
#' @param update_type deterministic or stochastic
#'
#' @examples
#' ## A single compartment and rate:
#' I <- 10; gamma <- 0.1; d_time <- 1
#' apply_rates(I, gamma, d_time, "deterministic")
#' ### Is the same as:
#' I * (1 - exp(-gamma*d_time))
#' ### Stochastic version uses rbinom:
#' apply_rates(I, gamma, d_time, "stochastic")
#'
#' ## A single compartment and vector of rates:
#' cull <- 0.01
#' apply_rates(I, c(gamma,cull), d_time, "deterministic")
#' ## Due to competing rates, this is NOT the same as:
#' I * (1 - exp(-gamma*d_time))
#' I * (1 - exp(-cull*d_time))
#' ### Stochastic version uses rmultinom:
#' apply_rates(I, c(gamma,cull), d_time, "stochastic")
#'
#' ## A vector of compartments and single rate:
#' E <- c(0,1,2); omega <- 0.1
#' apply_rates(E, omega, d_time, "deterministic")
#' ### Is the same as:
#' apply_rates(E[1], omega, d_time, "deterministic")
#' apply_rates(E[2], omega, d_time, "deterministic")
#' apply_rates(E[3], omega, d_time, "deterministic")
#' ### Is the same as:
#' E * (1 - exp(-omega*d_time))
#' ### Stochastic version uses rbinom:
#' apply_rates(E, omega, d_time, "stochastic")
#'
#' ## A vector of compartments and vector of rates:
#' repl <- 0.001
#' apply_rates(E, c(omega,repl), d_time, "deterministic")
#' apply_rates(E, c(omega,repl), d_time, "stochastic")
#'
#' ## If the rate vector has names, then these are carried over to the
#' ## output column names:
#' apply_rates(I, c(toR=gamma,toS=cull), d_time, "deterministic")
#'
#' ## If you only have a single rate (and/or single compartment) then you
#' ## probably want to remove the dimensions attribute before using, i.e.:
#' apply_rates(E, omega, d_time, "deterministic")[,1]
#' apply_rates(I, c(gamma,cull), d_time, "deterministic")[1,]
#' ## Or:
#' as.numeric(apply_rates(E, omega, d_time, "deterministic"))
#' as.numeric(apply_rates(I, c(gamma,cull), d_time, "deterministic"))
#'
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
      rmultinom(1, x, prob)[-1L,,drop=FALSE]
    }, numeric(length(prop))) ->
      rv
    if(length(compartments)==1L || length(prop)==1L){
      rv <- matrix(rv, nrow=length(compartments), ncol=length(prop))
    }else{
      rv <- t(rv)
    }
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
