#' Suggest a number of sub-compartments to obtain a suitable distribution of transition times
#'
#' @description
#' This is a utility function that takes a desired mean and standard deviation for
#' the distribution of transit times through a disease compartment, and suggests
#' a suitable rate and number of sub-compartments that matches the provided mean
#' and variance as closely as possible.
#'
#' @param mean the desired mean waiting time
#' @param sd the desired standard deviation in waiting time
#'
#' @returns the function returns a ggplot object illustrating the suggestions,
#' and prints the suggested number of sub-compartments and rates to screen
#'
#' @examples
#' # An average of 5 days, with sd of 2.5 days:
#' suggest_boxes(mean=5, sd=3)
#'
#' # The number of sub-compartments is invariant to transformation of the rate
#' # on different time scales, i.e. we get the same if specifying in hours:
#' suggest_boxes(mean=5*24, sd=3*24)
#'
#' @importFrom stats dgamma qgamma
#'
#' @export
suggest_boxes <- function(mean, sd){

  qassert(mean, "N1(0,)")
  qassert(sd, "N1(0,)")

  rate <- signif(1/mean, 2)

  #variance <- sd^2
  #shape <- mean^2/variance
  shape <- (mean/sd)^2

  options <- seq(floor(shape), ceiling(shape), by=1L)
  opt_shapes <- unique(options[options>0])
  vars <- opt_shapes/(opt_shapes*rate)^2

  yields <- paste0(opt_shapes, " subcompartments with overall rate ", rate, ", which yields:\n\tmean = ", round(1/rate, 2), "  (desired: ", round(mean,2), ")\n\tsd = ", round(sqrt(vars), 2), "  (desired: ", round(sd,2), ")\n", sep="")
  yields

  if(length(opt_shapes)==1L){
    cat("The closest match is ", yields, sep="")
    #cat("The closest match is Erlang(shape=", shapes, ", scale=", scales, "), which yields:\n\tmean = ", round(shapes/rates, 2), "  (desired: ", round(mean,2), ")\n\tsd = ", round(shapes^0.5/scales, 2), "  (desired: ", round(sd,2), ")\n", sep="")
  }else{
    cat("Possible matches are:\n", paste0("\t", str_replace_all(yields, "\n\t","\n\t\t")))
  }

  seq_along(opt_shapes) |>
    lapply(\(x){
      shape <- opt_shapes[x]
      scale <- mean/shape
      range <- qgamma(c(0.001,1-0.001), shape, scale=scale)
      tibble(
        `#Boxes: `=as.character(shape), Time=seq(range[1],range[2],length.out=1e3),
        `rate: `=as.character(rate)
      ) |>
        mutate(Density = dgamma(Time,shape,scale=scale))
    }) |>
    bind_rows() |>
    ggplot(aes(x=Time, y=Density, col=`#Boxes: `, lty=`rate: `)) +
    geom_line() +
    theme_light() +
    theme(legend.position="bottom")

}
