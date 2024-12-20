#' Discrete-time deterministic and stochastic SEIR models for multiple groups
#' @name multi_models
#'
#' @description
#' These wrapper functions create a multi-group model, run the model for the
#' indicated time steps, and then return a data frame representing the model
#' output. Note that starting values for S/E/I/R can be continuous or discrete
#' for the deterministic models, but must be discrete numbers for the stochastic
#' models.
#'
#' @param n_groups the number of groups to generate
#' @param beta_matrix a matrix representing between-group spread (must have nrow and ncol equal to n_groups)
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
#' @param transmission_within the within-group transmission type (frequency or density)
#' @param transmission_between the between-group transmission type (frequency or density)
#' @param iterations the number of iterations to run (stochastic models only)
#' @param d_time the desired time step (delta time)
#' @param max_time the desired maximum time point (must be greater than the time step)
#'
#' @examples
#' output <- multi_seir_det(2, matrix(0,nrow=2,ncol=2))
#' ggplot2::autoplot(output)
#'
#' output <- multi_seir_stoc(2, matrix(0,nrow=2,ncol=2), iterations=3L)
#' ggplot2::autoplot(output)
NULL


#' @rdname multi_models
#' @export
multi_seir_det <- function(n_groups, beta_matrix, S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_within=c("frequency","density"), transmission_between=c("density","frequency"), d_time=1, max_time=100){

  transmission_within <- match.arg(transmission_within)
  transmission_between <- match.arg(transmission_between)
  output <- multi_wrapper("deterministic", n_groups=n_groups, beta_matrix=beta_matrix, S=S, E=E, I=I, R=R, numE=numE, beta=beta, omega=omega, gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull, transmission_within=transmission_within, transmission_between=transmission_between, d_time=d_time, max_time=max_time, iterations=1L)

  class(output) <- c("ipdmr_dm", class(output))
  attr(output, "plot_caption") <- str_c("deterministic; ", n_groups, " groups")

  output
}

#' @rdname multi_models
#' @export
multi_seir_stoc <- function(n_groups, beta_matrix, S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_within=c("frequency","density"), transmission_between=c("density","frequency"), d_time=1, max_time=100, iterations=1L){

  transmission_within <- match.arg(transmission_within)
  transmission_between <- match.arg(transmission_between)
  output <- multi_wrapper("stochastic", n_groups=n_groups, beta_matrix=beta_matrix, S=S, E=E, I=I, R=R, numE=numE, beta=beta, omega=omega, gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull, transmission_within=transmission_within, transmission_between=transmission_between, d_time=d_time, max_time=max_time, iterations=iterations)

  class(output) <- c("ipdmr_sm", class(output))
  attr(output, "iterations") <- iterations
  attr(output, "plot_caption") <- str_c("stochastic; ", n_groups, " groups")

  output
}


multi_wrapper <- function(update_type, n_groups, beta_matrix, S=99, E=0, I=1, R=0, numE=3L, beta=0.25, omega=0.2, gamma=0.2, delta=0.05, vacc=0, repl=0, cull=0, transmission_within=c("frequency","density"), transmission_between=c("density","frequency"), d_time=1, max_time=100L, iterations=1L){

  ## Check arguments:
  transmission_within <- match.arg(transmission_within)
  transmission_between <- match.arg(transmission_between)
  qassert(n_groups, "X1(0,)")
  assert_matrix(beta_matrix, "numeric", any.missing=FALSE, nrows=n_groups, ncols=n_groups)
  qassert(iterations, "X1(0,)")
  qassert(d_time,"N1(0,)")
  qassert(max_time,"N1(0,)")

  model_pars <- list(
    S=S, E=E, I=I, R=R, numE=numE, beta=beta, omega=omega,
    gamma=gamma, delta=delta, vacc=vacc, repl=repl, cull=cull,
    transmission_type=transmission_within,
    group_number=seq_len(n_groups)
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
      md <- make_group(update_type=update_type, numE=x$numE, numI=1, numR=1, group_name=x[["group_number"]])
      x$numE <- NULL
      x$group_number <- NULL
      for(nm in names(x)){
        md[[nm]] <- x[[nm]]
      }
      md$save()
      md
    }) ->
    models

  ##  Set up the model:
  model <- BetweenGroupClass$new(models)
  model$beta_matrix <- beta_matrix
  model$transmission_between <- transmission_between
  model$save()

  ## Iterate:
  pblapply(seq_len(iterations), \(i){

    ## Reset and run:
    model$reset()
    model$run(max_time, d_time) |>
      as_tibble() |>
      mutate(Iteration = i) |>
      select("Iteration", everything())
  }) |>
    bind_rows() |>
    {\(x) if(update_type=="deterministic") select(x, -"Iteration") else x}()
}
