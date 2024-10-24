#' Title
#'
#' @param N
#' @param beta
#' @param type
#' @param init_I
#' @param time_points
#'
#' @importFrom checkmate qassert assert_number
#' @importFrom deSolve ode
#' @import stringr
#' @import dplyr
#'
#' @examples
#' si_continuous(N=10, type="density") |> ggplot2::autoplot()
#' si_continuous(N=10, type="frequency") |> ggplot2::autoplot()
#'
#' @export
si_continuous <- function(S=9, I=1, beta=0.05, type=c("frequency","density"), time_points=seq(0,21,by=0.1)){

  qassert(S, "N1[0,)")
  qassert(I, "N1[0,)")
  qassert(S+I, "N1(0,)")
  N <- S+I
  qassert(beta, "N1(0,)")
  qassert(time_points, "N+[0,)")

  type <- match.arg(type)

  model_dens <- function(times, y, parameters)
  {
    new_I <- parameters[["beta"]]*y["S"]*y["I"]
    list(c(-new_I, +new_I))
  }
  model_freq <- function(times, y, parameters)
  {
    new_I <- parameters[["beta"]]*y["S"]*y["I"]/parameters[["N"]]
    list(c(-new_I, new_I))
  }

  parameters <- list(N=N, beta=beta)

  ode(
    y = c("S"=S, "I"=I),
    times = time_points,
    func = if(type=="density") model_dens else model_freq,
    parms = parameters
  ) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(S = pmax(.data$S, 0), I = N-.data$S) |>
    select(Time=.data$time, everything()) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("continuous; ", type)
  return(output)
}


