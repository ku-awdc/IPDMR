#' Title
#'
#' @param S
#' @param I
#' @param R
#' @param beta
#' @param gamma
#' @param delta
#' @param transmission_type
#' @param time_step
#' @param max_time
#'
#' @export
sirs_det <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, transmission_type="frequency", time_step=1L, max_time=100L){

  model <- WithinGroupModel$new(model_type="sirs", update_type="deterministic", transmission_type=transmission_type, time_step=time_step)

  model$S <- S
  model$I <- I
  model$R <- R

  model$beta <- beta
  model$gamma <- gamma
  model$delta <- delta

  model$run(ceiling(max_time/time_step))

}


#' Title
#'
#' @param S
#' @param I
#' @param R
#' @param beta
#' @param gamma
#' @param delta
#' @param transmission_type
#' @param time_step
#' @param max_time
#'
#' @importFrom pbapply pblapply
#'
#' @export
sirs_stoc <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, iterations=1, transmission_type="frequency", time_step=1L, max_time=100L){

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="sirs", update_type="stochastic", transmission_type=transmission_type, time_step=time_step)

    model$S <- S
    model$I <- I
    model$R <- R

    model$beta <- beta
    model$gamma <- gamma
    model$delta <- delta

    model$run(ceiling(max_time/time_step)) |>
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

