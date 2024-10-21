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
sirs_stoch <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, iterations=1, transmission_type="frequency", time_step=1L, max_time=100L){

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="sirs", update_type="stochastic", transmission_type=transmission_type, time_step=time_step)

    model$S <- S
    model$I <- I
    model$R <- R

    model$beta <- beta
    model$gamma <- gamma
    model$delta <- delta

    model$run(ceiling(max_time/time_step)) |>
      mutate(Iteration = i) |>
      select(.data$Iteration, everything())

  }) |>
    bind_rows() ->
    output

  if(iterations==1) output$Iteration <- NULL
  attr(output, "iterations") <- iterations
  output
}

