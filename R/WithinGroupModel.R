`$` <- function(lhs,rhs){
  name <- deparse(substitute(rhs))
  if(is.environment(lhs) && !name %in% names(lhs)) stop(name, " not present in environment")
  do.call(base::`$`, args=list(lhs, name))
}
mlist <- function(..., .dummy=NULL){
  args <- list(...)
  if(any(names(args)=="")){
    an <- names(args)
    an[names(args)==""] <- "MISSING"
    stop("Empty names in member/method list: ", paste0(an,collapse=", "))
  }
  if(any(table(names(args))>1)){
    tn <- table(names(args))
    stop("Duplicated names in member/method list: ", paste0(names(tn)[tn>1], collapse=", "))
  }
  do.call(list, args)
}

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
      model_type = c("sirs"),
      update_type = c("deterministic","stochastic"),
      transmission_type = c("frequency","density"),
      time_step = 1L
    ){
      model_type <- match.arg(model_type)
      update_type <- match.arg(update_type)
      transmission_type <- match.arg(transmission_type)
      wgm_initialise(super, self, private, model_type, update_type, transmission_type, time_step)
    },

    update = function(
      time_step = self$time_step
    ){
      wgm_update_models[[private$.model_type]](super, self, private, time_step)
      invisible(self)
    },

    check_state = function(){

    },

    run = function(n_steps, time_step=self$time_step, include_current=self$time==0){
      wgm_run(super, self, private, n_steps, time_step, include_current)
    },

    .dummy=NULL
  ),

  private = mlist(
    .model_type = character(),
    .update_type = character(),
    .transmission_type = character(),
    .time_step = 1,

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

    time_step = function(value){
      if(missing(value)) return(private$.time_step)
      private$.time_step <- value
    },

    time = function() private$.time,

    state = function(){
      tibble(Time = self$time, S = self$S, I = self$I, R = self$R)
    },

    .dummy=function() NULL
  ),

  lock_class = TRUE

)

wgm_initialise <- function(super, self, private, model_type="hi", update_type, transmission_type, time_step){

  private$.model_type <- model_type
  private$.update_type <- update_type
  private$.transmission_type <- transmission_type
  private$.reset_N()
  self$time_step <- time_step

}

wgm_update_models <- list(
  sirs = function(super, self, private, time_step){

    self$check_state()

    S <- private$.S
    I <- private$.I
    R <- private$.R

    infctn <- if(private$.transmission_type=="frequency"){
      I/private$.N * self$time_step
    }else{
      I * self$time_step
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

    new_I <- cfun(S, (1 - exp(-private$.beta * infctn)))
    new_R <- cfun(I, (1 - exp(-private$.gamma)))
    new_S <- cfun(R, (1 - exp(-private$.delta)))

    private$.S <- S + new_S - new_I
    private$.I <- I + new_I - new_R
    private$.R <- R + new_R - new_S

    private$.time <- private$.time + time_step

    self$check_state()

  }
)

wgm_run <- function(super, self, private, n_steps, time_step, include_current){
  c(
    if(include_current) list(self$state),
    lapply(seq_len(n_steps), \(x){
      model$update(time_step)
      model$state
    })
  ) |>
    bind_rows() ->
    out
  cn <- if(private$.update_type=="deterministic") "ipdmr_dt" else "ipdmr_st"
  class(out) <- c(cn, class(out))
  out
}

# model <- WithinGroupModel$new(update_type="stochastic"); model$run(100) |> autoplot()
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
