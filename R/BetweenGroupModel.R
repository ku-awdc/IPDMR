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

    models = list(),
    beta_matrix = matrix(),

    update = function(){

      stopifnot(dim(self$beta_matrix)==length(self$models))
      diag(self$beta_matrix) <- 0.0

      trans_b <- colSums(self$beta_matrix * sapply(self$models, \(x) x$I))
      for(i in seq_along(self$models)){
        self$models[[i]]$trans_external <- trans_b[i]
        self$models[[i]]$update()
      }

      invisible(self)
    },

    run = function(n_steps, d_time=1){

      c(
        lapply(self$models, \(x) x$state),
        lapply(seq_len(n_steps), function(x){
          self$update()
          lapply(self$models, \(x) x$state) |> bind_rows()
        })
      ) |>
        bind_rows()
    },

    .dummy=NULL
  ),

  lock_class = TRUE
)

