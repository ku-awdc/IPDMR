if(length(find("mlists"))==0L) source("R/R6utils.R")

## TODO: call_method(mthd_name) function that takes a function mthd_name, extracts the args() as a string, gets rid of super, self, private, then adds a body to call substitute(mthd_name) with the args() intact (but no defaults i.e. remove anything after =)
## so we would have e.g. initialize = call_method(wgm_initialize)
## would this break roxygen docs?  I don't think so...


#' @title WithinGroupModel
#'
#' @description
#' General within-group model class
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
WithinGroupModel <- R6::R6Class(

  "WithinGroupModel",

  public = mlist(

    initialize = function(
      model_type = c("sir"),
      update_type = c("deterministic","stochastic"),
      transmission_type = c("frequency","density"),
      d_time = 1L
    ){
      model_type <- match.arg(model_type)
      update_type <- match.arg(update_type)
      transmission_type <- match.arg(transmission_type)
      wgm_initialise(super, self, private, model_type, update_type, transmission_type, d_time)
    },

    update = function(
    d_time = self$d_time
    ){
      wgm_update_models[[private$.model_type]](super, self, private, d_time)
      invisible(self)
    },

    check_state = function(){

    },

    run = function(n_steps, d_time=self$d_time, include_current=self$time==0){
      wgm_run(super, self, private, n_steps, d_time, include_current)
    },

    trans_external = 0.0,
    id = 0,

    .dummy=NULL
  ),

  private = mlist(
    .model_type = character(),
    .update_type = character(),
    .transmission_type = character(),
    .d_time = 1,

    .S = 99,
    .I = 1,
    .R = 0,
    .N = 100,
    .time = 0,
    .beta = 0.25,
    .gamma = 0.2,
    .delta = 0.05,

    .reset_N = function(){
      private$.N <- private$.S + private$.I + private$.R
    },

    .dummy = NULL
  ),

  active = list(
    S = function(value){
      if(missing(value)) return(private$.S)
      if(private$.update_type=="stochastic") stopifnot(value%%1==0)
      private$.S <- value
      private$.reset_N()
    },

    I = function(value){
      if(missing(value)) return(private$.I)
      if(private$.update_type=="stochastic") stopifnot(value%%1==0)
      private$.I <- value
      private$.reset_N()
    },

    R = function(value){
      if(missing(value)) return(private$.R)
      if(private$.update_type=="stochastic") stopifnot(value%%1==0)
      private$.R <- value
      private$.reset_N()
    },

    beta = function(value){
      if(missing(value)) return(private$.beta)
      stopifnot(value >= 0.0)
      private$.beta <- value
    },

    gamma = function(value){
      if(missing(value)) return(private$.gamma)
      stopifnot(value >= 0.0)
      private$.gamma <- value
    },

    delta = function(value){
      if(missing(value)) return(private$.delta)
      stopifnot(value >= 0.0)
      private$.delta <- value
    },

    d_time = function(value){
      if(missing(value)) return(private$.d_time)
      private$.d_time <- value
    },

    time = function() private$.time,

    state = function(){
      if(self$id==0){
        tibble(Time = self$time, S = self$S, I = self$I, R = self$R)
      }else{
        tibble(Model = self$id, Time = self$time, S = self$S, I = self$I, R = self$R)
      }
    },

    .dummy=function() NULL
  ),

  lock_class = TRUE

)

wgm_initialise <- function(super, self, private, model_type="hi", update_type, transmission_type, d_time){

  private$.model_type <- model_type
  private$.update_type <- update_type
  private$.transmission_type <- transmission_type
  private$.reset_N()
  self$d_time <- d_time

}

wgm_update_models <- list(
  sir = function(super, self, private, d_time){

    self$check_state()

    S <- private$.S
    I <- private$.I
    R <- private$.R

    infctn <- if(private$.transmission_type=="frequency"){
      I/private$.N
    }else{
      I
    }

    cfun <- if(private$.update_type=="stochastic"){
      function(size, prob){
        rbinom(1, size, prob)
      }
    }else{
      function(size, prob){
        size*prob
      }
    }

    new_I <- cfun(S, (1 - exp(-private$.beta * infctn * self$d_time - self$trans_external)))
    new_R <- cfun(I, (1 - exp(-private$.gamma * self$d_time)))
    new_S <- cfun(R, (1 - exp(-private$.delta * self$d_time)))

    private$.S <- S + new_S - new_I
    private$.I <- I + new_I - new_R
    private$.R <- R + new_R - new_S

    private$.time <- private$.time + d_time

    self$check_state()

  }
)

wgm_run <- function(super, self, private, n_steps, d_time, include_current){
  c(
    if(include_current) list(self$state),
    lapply(seq_len(n_steps), \(x){
      self$update(d_time)
      self$state
    })
  ) |>
    bind_rows() ->
    out

  if(private$.update_type=="deterministic"){

    class(out) <- c("ipdmr_dt", class(out))
    attr(out, "plot_caption") <- str_c("deterministic; discrete; ", private$.transmission_type)

  }else{

    out |>
      mutate(Iteration = 1L) |>
      select(Iteration, everything()) ->
      out

    class(out) <- c("ipdmr_st", class(out))
    attr(out,'iterations') <- 1L
    attr(out, "plot_caption") <- str_c("stochastic; discrete; ", private$.transmission_type)

  }

  out
}

# model <- WithinGroupModel$new(update_type="stochastic"); model$beta <- 0.5; model$run(100) |> autoplot()
#
# c(
#   list(model$state),
#   lapply(1:360, \(x){
#     model$update()
#     model$state
#   })
# ) |>
#   bind_rows() ->
#   out
# head(out)
#
#
