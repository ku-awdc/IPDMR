#' A simple continuous-time SI model
#'
#' @param S the starting number of susceptible animals (continuous number)
#' @param I the starting number of infected animals (continuous number)
#' @param beta the transmission rate per unit time (must be positive)
#' @param transmission_type either frequency or density
#' @param time_points a vector of time points to include in the data frame returned
#'
#' @importFrom checkmate qassert assert_number assert_matrix
#' @importFrom deSolve ode
#' @import stringr
#' @import dplyr
#'
#' @examples
#' si_continuous(S=9, I=1, transmission_type="density") |> ggplot2::autoplot()
#' si_continuous(S=9, I=1, transmission_type="frequency") |> ggplot2::autoplot()
#'
#' @export
si_continuous <- function(S=9, I=1, beta=0.05, transmission_type=c("frequency","density"), time_points=seq(0,21,by=0.1)){

  qassert(S, "N1[0,)")
  qassert(I, "N1[0,)")
  qassert(S+I, "N1(0,)")
  N <- S+I
  qassert(beta, "N1(0,)")
  qassert(time_points, "N+[0,)")

  type <- match.arg(transmission_type)

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

  class(output) <- c("ipdmr_ct", class(output))
  attr(output, "plot_caption") <- str_c("continuous; ", type)
  return(output)
}



#' A super-simple continuous-time two-compartment model
#'
#' @param rate the rate of progression from A to B per unit time (must be positive)
#' @param time_points a vector of time points to include in the data frame returned
#' @export
ab_continuous <- function(rate=0.1, time_points=seq(0,10,by=0.1)){

  model_fun <- function(times, y, parameters)
  {
    new_B <- parameters[["rate"]]*y["A"]
    list(c(-new_B, +new_B))
  }

  parameters <- list(rate=rate)

  ode(
    y = c("A"=1, "B"=0),
    times = time_points,
    func = model_fun,
    parms = parameters
  ) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(A = pmax(.data$A, 0)) |>
    select(Time=.data$time, A=.data$A) ->
    output

  class(output) <- c("ipdmr_ct", class(output))
  attr(output, "plot_caption") <- str_c("continuous time")
  return(output)
}


