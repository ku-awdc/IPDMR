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
seir_det <- function(S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_type="frequency", d_time=1L, max_time=100L){

  model <- WithinGroupModel$new(model_type="seir", update_type="deterministic", transmission_type=transmission_type, numE=numE, d_time=d_time)

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
seir_stoc <- function(S=99, E=0, I=1, R=0, num_E=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0,iterations=1, transmission_type="frequency", d_time=1L, max_time=100L){

  pblapply(seq_len(iterations), \(i){

    model <- WithinGroupModel$new(model_type="seir", update_type="stochastic", transmission_type=transmission_type, num_E=num_E, d_time=d_time)

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


#' @export
multi_seir_det <- function(n_groups, beta_matrix, S=99, E=0, I=1, R=0, num_E=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, d_time=1, max_time=100L){

  output <- multi_wrapper("deterministic", n_groups=n_groups, beta_matrix=beta_matrix, S=S, E=E, I=I, R=R, num_E=num_E, beta=beta, omega=omega, gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull, iterations=1, d_time=d_time, max_time=max_time)

  class(output) <- c("ipdmr_dm", class(output))
  attr(output, "plot_caption") <- str_c("deterministic; ", n_groups, " groups")

  output
}


#' @export
multi_seir_stoc <- function(n_groups, beta_matrix, S=99, E=0, I=1, R=0, num_E=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, iterations=1, d_time=1, max_time=100L){

  output <- multi_wrapper("stochastic", n_groups=n_groups, beta_matrix=beta_matrix, S=S, E=E, I=I, R=R, num_E=num_E, beta=beta, omega=omega, gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull, iterations=iterations, d_time=d_time, max_time=max_time)

  class(output) <- c("ipdmr_sm", class(output))
  attr(output, "iterations") <- iterations
  attr(output, "plot_caption") <- str_c("stochastic; ", n_groups, " groups")

  output
}


multi_wrapper <- function(update_type, n_groups, beta_matrix, S=99, E=0, I=1, R=0, num_E=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, iterations=1, d_time=1, max_time=100L){

  ## Check arguments:
  qassert(n_groups, "X1(0,)")
  assert_matrix(beta_matrix, "numeric", any.missing=FALSE, nrows=n_groups, ncols=n_groups)
  qassert(iterations, "X1(0,)")
  qassert(d_time,"N1(0,)")
  qassert(max_time,"N1(0,)")

  model_pars <- list(
    S=S, E=E, I=I, R=R, num_E=num_E, beta=beta, omega=omega,
    gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull
  )

  ## Check parameters then build within-group models:
  names(model_pars) |>
    set_names() |>
    lapply(\(x){
      if(length(model_pars[[x]])==1L) model_pars[[x]] <- rep(model_pars[[x]], n_groups)
      if(!length(model_pars[[x]])==n_groups) stop("Invalid length for ", x, ": ", length(model_pars[[x]]), " (should be 1 or ", n_groups, ")", call.=FALSE)
      model_pars[[x]]
    }) |>
    as_tibble() |>
    rowwise() |>
    group_split() |>
    set_names(str_c("Group_", seq_len(n_groups) |> format() |> str_replace_all(" ", "0"))) |>
    lapply(\(x){
      md <- WithinGroupModel$new("seir", update_type, "frequency", num_E=x$num_E, d_time=d_time)
      x$num_E <- NULL
      for(nm in names(x)){
        md[[nm]] <- x[[nm]]
      }
      md
    }) ->
    models

  ## Set up between-group model (we do a deep clone later, so no need to clone here):
  bgmaster <- BetweenGroupModel$new(models, d_time=d_time, clone_models=FALSE)
  bgmaster$beta_matrix <- beta_matrix

  ## Iterate:
  pblapply(seq_len(iterations), \(i){

    model <- bgmaster$clone(deep=TRUE)
    model$run(ceiling(max_time/d_time)) |>
      as_tibble() |>
      mutate(Iteration = i) |>
      select("Iteration", everything())
  }) |>
    bind_rows() |>
    {\(x) if(update_type=="deterministic") select(x, -.data$Iteration) else x}()
}
