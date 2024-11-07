#' Suggest a number of sub-compartments to obtain a suitable distribution of transition times
#'
#' @description
#' This is a utility function that takes a desired mean and standard deviation for
#' the distribution of transit times through a disease compartment, and suggests
#' a suitable rate and number of sub-compartments that matches the provided mean
#' and variance as closely as possible.  If a d_time argument is provided, then
#' the illustrative plot of the Erlang (gamma) distribution will be overlaid
#' with a negative binomial.
#'
#' @param mean
#' @param sd
#' @param d_time
#'
#' @returns the function returns a ggplot object illustrating the suggestions,
#' and prints the suggested number of sub-compartments and rates to screen
#'
#' @examples
#' suggest_boxes(mean=5, sd=1)
#' suggest_boxes(mean=5, sd=3, d_time=0.1)
#'
#' @export
suggest_boxes <- function(mean, sd, d_time){

  if(missing(d_time)) d_time <- NA_real_

  stop("IMPLEMENT ME")

  variance <- sd^2
  shape <- mean^2/variance^2

  ## TODO: calculate 99% CI from qgamma as from/to
  curve(dgamma(x, shape, scale=mean/shape), from=max(0,mean-5*sd), to=mean+5*sd)

  options <- c(floor(shape), ceiling(shape), ceiling(shape)+1)
  shapes <- unique(options[options>0])
  rates <- round(1/mean, 2)
  scales <- shapes/mean

  yields <- paste0(shapes, " subcompartments with overall rate ", rates, ", which yields:\n\tmean = ", round(1/rates, 2), "  (desired: ", round(mean,2), ")\n\tsd = ", round(shapes^0.5*scales, 2), "  (desired: ", round(sd,2), ")\n", sep="")
  yields

  if(length(shapes)==1L){
    cat("The closest match is ", yields, sep="")
    #cat("The closest match is Erlang(shape=", shapes, ", scale=", scales, "), which yields:\n\tmean = ", round(shapes/rates, 2), "  (desired: ", round(mean,2), ")\n\tsd = ", round(shapes^0.5/scales, 2), "  (desired: ", round(sd,2), ")\n", sep="")
  }else{
    cat("Possible matches are:\n", paste0("\t", str_replace_all(yields, "\n\t","\n\t\t")))
  }

  stop("Add mandatory N property to within-gp and show/add freq-based for beta matrix, and allow dens-based for multi wrapper")

}
