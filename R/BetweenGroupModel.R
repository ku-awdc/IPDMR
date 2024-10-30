if(length(find("mlists"))==0L) source("R/R6utils.R")

#' @title BetweenGroupModel
#'
#' @description
#' General between-group model class
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
BetweenGroupModel <- R6::R6Class(

  "BetweenGroupModel",

  public = mlist(

    initialize = function(models, d_time=0.1, clone_models=FALSE){

      if(inherits(models, "R6")){
        models <- list(models)
      }
      if(!is.list(models) || length(list)==0L || any(!sapply(models, \(x) inherits(x, "R6")))){
        stop("The input models must be a list of 1 or more R6 objects")
      }
      if(!all(sapply(models, \(x){
        all(c("I","state","trans_external","update") %in% names(x)) &&
          length(formals(x$update))==1L
      }))){
        stop("All input models must have public fields (or active bindings) for I, state and trans_external, as well as a public update method with a single argument (d_time)")
      }

      if(clone_models){
        private$.models <- apply(models, \(x) x$clone())
      }else{
        private$.models <- models
      }

      private$.ngroups <- length(models)

      qassert(d_time, "N1(0,)")
      private$.d_time <- d_time

      private$.beta_matrix <- matrix(0, ncol=private$.ngroups, nrow=private$.ngroups)
      private$.time <- 0

    },

    update = function(d_time=private$.d_time){

      qassert(d_time, "N1(0,)")
      stopifnot(dim(private$.beta_matrix)==private$.ngroups)

      trans_b <- colSums(private$.beta_matrix * sapply(private$.models, \(x) x$I))
      for(i in seq_along(private$.models)){
        private$.models[[i]]$trans_external <- trans_b[i]
        private$.models[[i]]$update()
      }

      invisible(self)
    },

    run = function(n_steps, d_time=private$.d_time, iterations=1L, include_current=self$time==0){

      c(
        if(include_current) lapply(private$.models, \(x) x$state) else list(),
        lapply(seq_len(n_steps), function(x){
          self$update()
          lapply(private$.models, \(x) x$state) |> bind_rows()
        })
      ) |>
        bind_rows() |>
        ## In case not all models have the same compartments:
        mutate(across(is.numeric, \(x) replace_na(x, 0))) ->
        out
      class(out) <- c("ipdmr_mt", class(out))
      out

    },

    .dummy=NULL
  ),

  private = mlist(

    .beta_matrix = matrix(),
    .ngroups = numeric(),
    .models = list(),
    .d_time = numeric(),
    .time = numeric(),

    .dummy=NULL
  ),

  active = mlist(

    beta_matrix = function(value){
      if(missing(value)) return(private$.beta_matrix)
      stopifnot(value >= 0.0, dim(value)==private$.ngroups)
      private$.beta_matrix <- value
    },

    models = function(){
      private$.models
    },

    time = function(){
      private$.time
    },

    .dummy=NULL
  ),

  lock_class = TRUE
)

