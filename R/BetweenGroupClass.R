if(length(find("mlists"))==0L) source("R/R6utils.R")

#' Encapsulated OO class for between-group spread models
#'
#' @description
#' This is a general-purpose between-group model spread class that takes a list
#' of within-group spread models as input, and provides methods to update/run
#' (and save/reset) these models while taking care of between-group
#' (density-dependent) spread by calculating and inputting trans_external for
#' each of the groups as appropriate based on the beta_matrix, which can be
#' modified as needed (see examples).
#'
#' The input models can either be created using the \link{make_group} function,
#' or they can be your own R6 (or C++) class that implement the necessary methods
#' and fields (i.e. reset, save, update, I, state, and trans_external). The
#' models need not all be the same subtype, i.e. you can mix SIR/SEIR/SI models
#' and even R6 and C++ implementations.
#'
#' Note that the between-group model stores a reference to (i.e. shallow copy of)
#' the input models, so you can see the state of the input groups being modified
#' as the between-group model runs. You can also change the state of these groups
#' at any time between updates/runs, either by modifying the groups used to create
#' the between-group model class or by using the (read-only) active binding groups
#' to retrieve one or more reference to a group. Of course you can also change
#' the beta_matrix at any time between updates/runs.
#'
#'
#' @examples
#' ## Create a list of groups:
#' groups <- list(
#'   make_group("stochastic", group_name="Group1"),
#'   make_group("stochastic", group_name="Group2")
#' )
#'
#' ## Modify the starting conditions as needed:
#' groups[[1]]$S <- 100
#' groups[[1]]$I <- 0
#' groups[[2]]$S <- 99
#' groups[[2]]$I <- 1
#' ## etc
#'
#' ## Create a between-group model:
#' model <- BetweenGroupClass$new(groups)
#'
#' ## Modify beta_matrix as needed:
#' model$beta_matrix <- matrix(c(0,0.1,0.01,0), nrow=2, ncol=2)
#'
#' ## Then run the model:
#' output <- model$run(21, d_time=0.1)
#' ggplot2::autoplot(output)
#'
#' ## Change the beta_matrix to e.g. remove transmission from Group1->Group2:
#' model$beta_matrix <- matrix(c(0,0,0.01,0), nrow=2, ncol=2)
#' ## And run for more time points:
#' output2 <- model$run(21, d_time=0.1)
#' ggplot2::autoplot(output2)
#'
#' ## To combine time points:
#' all_output <- dplyr::bind_rows(output, output2)
#' summary(all_output)
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#'
#' @export
BetweenGroupClass <- R6::R6Class(

  "BetweenGroupClass",

  public = mlist(

    #' @description
    #' Create a new between-group model
    #' @param groups a list of within-group models as either R6 or C++ classes (see e.g. \link{make_group})
    #' @return A new between-group model object
    initialize = function(groups){

      ## Allow a single model as non-list:
      if(!inherits(groups, "list")){
        groups <- list(groups)
      }
      if(!is.list(groups) || length(list)==0L || any(!sapply(groups, \(x) inherits(x, "R6")))){
        # TODO: fails for C++
        #stop("The input models must be a list of 1 or more R6 objects")
      }
      if(!all(sapply(groups, \(x){
        all(c("I","state","trans_external","update") %in% names(x)) &&
          length(formals(x$update))==1L
      }))){
        # TODO: fails for C++
        #stop("All input models must have public fields (or active bindings) for I, state and trans_external, as well as a public update method with a single argument (d_time)")
      }

      ## TODO: set value correctly
      private$.allcpp <- FALSE

      private$.groups <- groups
      private$.ngroups <- length(groups)

      private$.beta_matrix <- matrix(0, ncol=private$.ngroups, nrow=private$.ngroups)
      private$.time <- 0

      self$save()

    },

    #' @description
    #' Update the state of each group for a single time point, including
    #' setting trans_external based on the value of the beta_matrix and vector
    #' of I obtained from the groups (density-dependent spread).
    #' @param d_time the desired time step (delta time)
    #' @return self, invisibly
    update = function(d_time){

      qassert(d_time, "N1(0,)")
      stopifnot(dim(private$.beta_matrix)==private$.ngroups)

      ## TODO: C++ function if private$.allcpp

      trans_b <- colSums(private$.beta_matrix * sapply(private$.groups, \(x) x$I))
      for(i in seq_along(private$.groups)){
        private$.groups[[i]]$trans_external <- trans_b[i]
        private$.groups[[i]]$update(d_time)
      }

      invisible(self)
    },

    #' @description
    #' Instruct each group to save the current state and parameter values for later retrieval using reset()
    #' @return self, invisibly
    save = function(){
      lapply(private$.groups, \(x) x$save())
      private$.time <- 0
    },

    #' @description
    #' Instruct each group to reset the current state and parameter values to their last saved state
    #' @return self, invisibly
    reset = function(){
      lapply(private$.groups, \(x) x$reset())
      private$.time <- 0
    },

    #' @description
    #' Update the state of each group for a given number of time points
    #' @param add_time the additional time to add to the current time of the model
    #' @param d_time the desired time step (delta time)
    #' @return a data frame of the model state at each (new) time point
    run = function(add_time, d_time){

      qassert(add_time, "N1(0,)")
      qassert(d_time, "N1(0,)")
      stopifnot(dim(private$.beta_matrix)==private$.ngroups)

      ## TODO: C++ function if private$.allcpp

      c(
        if(self$time==0){
          list(self$state)
        }else{
          list()
        },
        lapply(seq(self$time+d_time, self$time+add_time, by=d_time), function(x){
          self$update(d_time)$state
        })
      ) |>
        bind_rows() ->
        out
      class(out) <- c("ipdmr_dm", class(out))
      attr(out, "ngroups") <- private$.ngroups
      out

    },

    #' @description
    #' Update the state of each group for several time points, until a stopping
    #' criterion is met.
    #' NOTE: this method is a stub: it has not yet been implemented!
    #' @param criterion_fun a function taking an input data frame of states and returning a logical scalar indicating if the simulation should be stopped (TRUE) or continue to be
    #' updated (FALSE)
    #' @param d_time the desired time step (delta time)
    #' @return a data frame of the model state at each (new) time point
    run_until = function(criterion_fun, d_time){
      stop("This method is not yet implemented")
      ## TODO: criterion should be a function accepting a data frame of states
      ## and returning a logical i.e. continue/stop
    },

    #' @description
    #' Print method showing the number of groups and current time
    #' @return self, invisibly
    print = function(){
      cat("A between-group model with ", private$.ngroups, " within-group models\n", sep="")
      cat("[Current time point = ", private$.time, "]", sep="")
      invisible(self)
    },

    .dummy=NULL
  ),

  private = mlist(

    .allcpp = FALSE,
    .beta_matrix = matrix(),
    .ngroups = numeric(),
    .groups = list(),
    .time = numeric(),

    .dummy=NULL
  ),

  active = mlist(

    #' @field beta_matrix a matrix representing between-group (density-dependent) spread
    beta_matrix = function(value){
      if(missing(value)) return(private$.beta_matrix)
      stopifnot(value >= 0.0, dim(value)==private$.ngroups)
      private$.beta_matrix <- value
    },

    #' @field groups a list of the internally-stored within-group models (this can be used to directly check and/or modify their state and/or parameters)
    groups = function(){
      private$.groups
    },

    #' @field time the current time point
    time = function(){
      private$.time
    },

    #' @field state a data frame of the current state of each group
    state = function(){
      lapply(private$.groups, function(x) x$state) |>
        bind_rows() |>
        ## In case not all models have the same compartments:
        mutate(across(is.numeric, \(x) replace_na(x, 0))) |>
        mutate(GroupIndex = str_c("Gp", format(seq_along(private$.groups)) |> str_replace_all(" ", "0"))) |>
        select("GroupIndex", everything())
    },

    .dummy=NULL
  ),

  lock_class = TRUE
)

