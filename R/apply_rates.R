#' Title
#'
#' @param compartments
#' @param rates
#' @param d_time
#' @param update_type deterministic (default) or stochastic
#'
#' @export
apply_rates <- function(compartments, rates, d_time, update_type=c("deterministic","stochastic")){

  qassert(compartments, "N+[0,)")
  qassert(rates, "N+[0,)")
  qassert(d_time, "N1(0,)")
  update_type <- match.arg(update_type)

  ## Generate proportions and remove overlap:
  leave <- 1-exp(-sum(rates))
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
