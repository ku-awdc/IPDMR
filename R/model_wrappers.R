#' Discrete-time deterministic and stochastic SIR and SEIR models
#' @name group_models
#'
#' @description
#' These wrapper functions create a single-group model, run the model for the
#' indicated time steps, and then return a data frame representing the model
#' output. Note that starting values for S/E/I/R can be continuous or discrete
#' for the deterministic models, but must be discrete numbers for the stochastic
#' models.
#'
#'
#' @param S the starting number of susceptible animals
#' @param E the total starting number of exposed/latent animals (this will be distributed among the sub-compartments of E when the model is initialised)
#' @param I the starting number of infected animals
#' @param R the starting number of recovered animals
#' @param numE the number of sub-compartments within the exposed/latent state
#' @param beta the transmission rate parameter per unit time (must be positive)
#' @param omega the latent progression rate parameter per unit time (must be positive)
#' @param gamma the recovery rate parameter per unit time (must be positive)
#' @param delta the reversion rate parameter per unit time (must be positive)
#' @param vacc the vaccination rate parameter per unit time (must be positive)
#' @param repl the replacement rate parameter per unit time (must be positive)
#' @param cull the targeted culling rate parameter per unit time (must be positive)
#' @param iterations the number of iterations to run (stochastic models only)
#' @param transmission_type either frequency or density
#' @param d_time the desired time step (delta time)
#' @param max_time the desired maximum time point (must be greater than the time step)
#'
#' @importFrom pbapply pblapply
NULL

#' @rdname group_models
#' @export
sir_det <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, transmission_type="frequency", d_time=1, max_time=100){

  qassert(transmission_type, "S1")

  model <- WithinGroupModel$new(model_type="sir", update_type="deterministic", transmission_type=transmission_type, d_time=d_time)

  model$S <- S
  model$I <- I
  model$R <- R

  model$beta <- beta
  model$gamma <- gamma
  model$delta <- delta

  model$run(ceiling(max_time/d_time))

}


#' @rdname group_models
#' @export
seir_det <- function(S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_type=c("frequency","density"), d_time=1, max_time=100){

  #qassert(transmission_type, "S1")
  transmission_type <- match.arg(transmission_type)

  model <- make_group(update_type="deterministic", numE=numE, group_name=NA_character_, model_type="SEIR", implementation="C++")
  model$transmission_type <- transmission_type

  if(FALSE){
    model <- WithinGroupModel$new(model_type="seir", update_type="deterministic", transmission_type=transmission_type, numE=numE, d_time=d_time)
  }

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

  model$run(ceiling(max_time/d_time), d_time) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("deterministic; discrete; ", transmission_type)

  output
}


#' @rdname group_models
#' @export
sir_stoc <- function(S=99, I=1, R=0, beta=0.25, gamma=0.2, delta=0.05, iterations=1, transmission_type="frequency", d_time=1, max_time=100){

  qassert(transmission_type, "S1")

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="sir", update_type="stochastic", transmission_type=transmission_type, d_time=d_time)

    model$S <- S
    model$I <- I
    model$R <- R

    model$beta <- beta
    model$gamma <- gamma
    model$delta <- delta

    model$run(ceiling(max_time/d_time), d_time) |>
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


#' @rdname group_models
#' @export
seir_stoc <- function(S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0,iterations=1, transmission_type="frequency", d_time=1, max_time=100){

  #qassert(transmission_type, "S1")
  transmission_type <- match.arg(transmission_type)

  model <- make_group(update_type="stochastic", numE=numE, group_name=NA_character_, model_type="SEIR", implementation="C++")
  model$transmission_type <- transmission_type
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
  model$save()

  pblapply(seq_len(iterations), \(i){

    #model <- WithinGroupModel$new(model_type="seir", update_type="stochastic", transmission_type=transmission_type, num_E=num_E, d_time=d_time)

    model$reset()
    model$run(ceiling(max_time/d_time), d_time) |>
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

