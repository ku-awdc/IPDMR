#' A simple discrete-time SI model
#'
#' @param S the starting number of susceptible animals (continuous number)
#' @param I the starting number of infected animals (continuous number)
#' @param beta the transmission rate per unit time (must be positive)
#' @param transmission_type either frequency or density
#' @param d_time the desired time step (delta time)
#' @param max_time the desired maximum time point (must be greater than the time step)
#'
#' @importFrom checkmate qassert assert_number
#'
#' @examples
#' si_discrete(S=9, I=1, transmission_type="density") |> ggplot2::autoplot()
#' si_discrete(S=9, I=1, transmission_type="frequency") |> ggplot2::autoplot()
#'
#' @export
si_discrete <- function(S=9, I=1, beta=0.05, transmission_type=c("frequency","density"), d_time=1/24, max_time=21){

  qassert(S, "N1[0,)")
  qassert(I, "N1[0,)")
  qassert(S+I, "N1(0,)")
  N <- S+I
  qassert(beta, "N1(0,)")
  qassert(d_time, "N1[0,)")
  qassert(max_time, "N1[0,)")
  assert_number(max_time, lower=d_time)

  ntime <- ceiling(max_time/d_time)+1L
  type <- match.arg(transmission_type)

  model <- make_group("deterministic", numE=0L, numI=1L, numR=0L)
  model$transmission_type <- transmission_type
  model$S <- S
  model$I <- I
  model$beta <- beta
  model$gamma <- 0
  model$repl <- 0
  model$cull <- 0
  output <- model$run(max_time, d_time)
  rm(model)

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("discrete; ", type, "; d_time=", round(d_time,3))
  return(output)

  ### OLDER CODE BELOW HERE

  i_dens <- function(i) 1 - exp(-beta * i)
  i_freq <- function(i) 1 - exp(-beta * i/N)
  i_fun <- if(type=="density") i_dens else i_freq

  output <- matrix(NA_real_, nrow=ntime+1L, ncol=3L, dimnames=list(NULL, c("Time","S","I")))
  time <- 0
  output[1L,"Time"] <- time
  output[1L,"S"] <- S
  output[1L,"I"] <- I
  for(row in seq_len(ntime)){
    time <- time + d_time
    newI <- S*i_fun(I) * d_time
    S <- S - newI
    I <- I + newI
    output[row+1L,"Time"] <- time
    output[row+1L,"S"] <- S
    output[row+1L,"I"] <- I
  }

  output |>
    as_tibble() |>
    select("Time", everything()) |>
    filter(.data$Time <= max_time) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("discrete; ", type, "; d_time=", round(d_time,3))
  return(output)
}



#' A super-simple discrete-time two-compartment model
#'
#' @param rate the rate of progression from A to B per unit time (must be positive)
#' @param d_time the desired time step (delta time)
#' @param max_time the desired maximum time point (must be greater than the time step)
#' @export
ab_discrete <- function(rate=0.1, d_time=0.1, max_time=10){

  qassert(rate, "N1(0,)")
  qassert(d_time, "N1[0,)")
  qassert(max_time, "N1[0,)")
  assert_number(max_time, lower=d_time)

  ntime <- ceiling(max_time/d_time)+1L
  A <- 1

  output <- matrix(NA_real_, nrow=(ntime*2)+1L, ncol=2L, dimnames=list(NULL, c("Time","A")))
  time <- 0

  if(FALSE){
    A <- A - (A * (1 - exp(-rate*d_time)))
  }

  output[1,"Time"] <- 0
  output[1,"A"] <- A

  for(row in seq(2,ntime*2,by=2)){
    time <- time + d_time
    output[row,"Time"] <- time-1e-6
    output[row,"A"] <- A

    newA <- A * (1 - exp(-rate*d_time))
    A <- A - newA
    output[row+1,"Time"] <- time
    output[row+1,"A"] <- A
  }

  output |>
    as_tibble() |>
    select("Time", everything()) |>
    filter(.data$Time <= max_time) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("discrete; d_time=", round(d_time,3))
  return(output)
}


