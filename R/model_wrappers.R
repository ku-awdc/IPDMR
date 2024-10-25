#' Title
#'
#' @param S
#' @param I
#' @param R
#' @param beta
#' @param gamma
#' @param delta
#' @param transmission_type
#' @param d_time
#' @param max_time
#'
#' @importFrom pbapply pblapply
#'
#' @export
sir_det <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, transmission_type="frequency", d_time=1L, max_time=100L){

  model <- WithinGroupModel$new(model_type="sir", update_type="deterministic", transmission_type=transmission_type, d_time=d_time)

  model$S <- S
  model$I <- I
  model$R <- R

  model$beta <- beta
  model$gamma <- gamma
  model$delta <- delta

  model$run(ceiling(max_time/d_time))

}


#' @export
seir_det <- function(S=99, E=c(0,0,0), I=1, R=0, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_type="frequency", d_time=1L, max_time=100L){

  model <- WithinGroupModel$new(model_type="seir", update_type="deterministic", transmission_type=transmission_type, d_time=d_time)

  model$S <- S
  model$E <- E
  model$I <- I
  model$R <- R

  model$beta <- beta
  model$omega <- omega
  model$gamma <- gamma
  model$delta <- delta
  model$vacc <- vacc
  model$repl <- repl
  model$cull <- cull

  model$run(ceiling(max_time/d_time))
}


#' @export
sir_stoc <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, iterations=1, transmission_type="frequency", d_time=1L, max_time=100L){

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="sir", update_type="stochastic", transmission_type=transmission_type, d_time=d_time)

    model$S <- S
    model$I <- I
    model$R <- R

    model$beta <- beta
    model$gamma <- gamma
    model$delta <- delta

    model$run(ceiling(max_time/d_time)) |>
      as_tibble() |>
      mutate(Iteration = i) |>
      select("Iteration", everything())

  }) |>
    bind_rows() ->
    output

  class(output) <- c("ipdmr_st", class(output))
  attr(output, "iterations") <- iterations
  attr(output, "plot_caption") <- str_c("stochastic; discrete; ", transmission_type)

  output
}


#' @export
seir_stoc <- function(S=99, E=c(0,0,0), I=1, R=0, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0,iterations=1, transmission_type="frequency", d_time=1L, max_time=100L){

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="seir", update_type="stochastic", transmission_type=transmission_type, d_time=d_time)

    model$S <- S
    model$E <- E
    model$I <- I
    model$R <- R

    model$beta <- beta
    model$omega <- omega
    model$gamma <- gamma
    model$delta <- delta
    model$vacc <- vacc
    model$repl <- repl
    model$cull <- cull

    model$run(ceiling(max_time/d_time)) |>
      as_tibble() |>
      mutate(Iteration = i) |>
      select("Iteration", everything())
  }) |>
    bind_rows() ->
    output

  class(output) <- c("ipdmr_st", class(output))
  attr(output, "iterations") <- iterations
  attr(output, "plot_caption") <- str_c("stochastic; discrete; ", transmission_type)

  output
}
