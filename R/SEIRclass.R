#' Encapsulated OO classes for within-group SEIR models
#'
#' Both R6 and C++ implementations are available - see make_seir_model for
#' a wrapper function
#'
#' @import R6
#' @export
SEIRclass <- R6::R6Class("SEIRclass",

  public = mlist(

    #' @description
    #' Create a new within-group model
    #' @param update_type
    #' @param numE
    #' @param d_time
    #' @param group_name
    #' @return A new within-group model object
    initialize = function(update_type = c("deterministic","stochastic"), numE = 3L, d_time=0.1, group_name=NA_character_){
      ## Checking and partial matching:
      update_type <- match.arg(update_type)
      private$.update_type <- update_type

      ## Check specified number of E, then update E_ accordingly
      ## (it will be zero-initialised):
      qassert(numE, "X1(0,)")
      private$.numE <- as.integer(numE)
      private$.E <- numeric(private$.numE)
      private$reset_N()

      ## Check and save default d_time:
      qassert(d_time, "N1(0,)")
      private$.d_time <- d_time

      ## Save group name:
      qassert(group_name, "s1")
      private$.group_name <- group_name

      ## Set a default value for transmission type
      ## (this could also be an argument to initialize):
      self$transmission_type <- "frequency"

      ## Run save method:
      self$save()
      ## Run reset method to test that it works:
      self$reset()
    },

    ## New save method:
    save = function(){
      ## Shortcut to save all relevant public fields/active bindings:
      flds <- c("S","E","I","R","beta","omega","gamma","delta",
        "vacc","repl","cull","trans_external")
      private$.saved_public <- lapply(flds, function(n) self[[n]])
      ## Sanity check to make sure everything has been set to a non-null value:
      stopifnot(!sapply(private$.saved_public, is.null))
      ## Copy names:
      names(private$.saved_public) <- flds

      ## Shortcut to save all relevant private fields:
      flds <- c(".time")
      private$.saved_private <- lapply(flds, function(n) private[[n]])
      ## Sanity check to make sure everything has been set to a non-null value:
      stopifnot(!sapply(private$.saved_private, is.null))
      ## Copy names:
      names(private$.saved_private) <- flds

      invisible(self)
    },

    ## New reset method:
    reset = function(){
      ## Reset public fields and/or active bindings:
      for(n in names(private$.saved_public)){
        self[[n]] <- private$.saved_public[[n]]
      }
      ## Reset private fields:
      for(n in names(private$.saved_private)){
        private[[n]] <- private$.saved_private[[n]]
      }

      ## Re-run reset_N separately:
      private$reset_N()

      invisible(self)
    },

    ## Modified update method:
    update = function(d_time = private$.d_time){

      assert_number(d_time, lower=0)

      ## New safety feature:
      private$check_state()

      ## Note that transmission rate is now calculated by a separate private
      ## method, and that we pass the update_type directly to apply_rates:
      leave_S <- apply_rates(private$.S,
        c(private$get_transmission_rate(), private$.vacc),
        d_time,
        update_type = private$.update_type)
      leave_E <- apply_rates(private$.E,
        c(private$.omega*private$.numE, private$.repl),
        d_time,
        update_type = private$.update_type)
      leave_I <- apply_rates(private$.I,
        c(private$.gamma, private$.repl + private$.cull),
        d_time,
        update_type = private$.update_type)
      leave_R <- apply_rates(private$.R,
        private$.delta+private$.repl,
        d_time,
        update_type = private$.update_type)

      ## As for exercise 3B:
      private$.S <- private$.S + leave_R[,1] + sum(leave_E[,2]) + leave_I[,2] - sum(leave_S)
      private$.E <- private$.E + c(leave_S[,1], leave_E[-private$.numE,1]) - apply(leave_E,1,sum)
      private$.I <- private$.I + leave_E[private$.numE,1] - sum(leave_I)
      private$.R <- private$.R + leave_I[,1] + leave_S[,2] - sum(leave_R)

      private$.time <- private$.time + d_time

      ## New safety feature:
      private$check_state()

      invisible(self)
    },

    ## Exactly the same as before:
    run = function(add_time, d_time=private$.d_time, include_current=self$time==0){
      c(
        if(include_current) list(self$state) else list(),
        lapply(seq(self$time+d_time, self$time+add_time, by=d_time), function(x){
          self$update(d_time)$state
        })
      ) |>
        bind_rows()
    },

    get_state = function(){
      self$state
    },

    ## New compact (ish!) print method:
    print = function(){
      if(is.na(private$.group_name)){
        cat("An SEIR model with ")
      }else{
        cat("An SEIR model with identifier/name '", private$.group_name, "' and ", sep="")
      }
      cat("the following properties:\n\t",
        "S/E/I/R (N) = ", self$S, "/", self$E, "/", self$I, "/", self$R, " (", self$N, ")\n\t",
        "beta/omega/gamma/delta = ", self$beta, "/", self$omega, "/", self$gamma, "/", self$delta, "\n\t",
        "vacc/repl/cull = ", self$vacc, "/", self$repl, "/", self$cull, "\n\t",
        "E compartments = ", private$.numE, "\n\t",
        "external transmission = ", private$.trans_external, "\n\t",
        "update type = ", private$.update_type, "\n\t",
        "transmission type = ", private$.transmission_type, "\n\t",
        "current time = ", self$time, "\n\t",
        sep="")
    }
  ),

  private = mlist(

    ## Private fields - each has a trailing underscore to avoid name clashes:
    .S = 99,
    .E = numeric(), # The length of this is set at initialisation
    .I = 1,
    .R = 0,
    .N = numeric(), # Set by reset_N method
    .update_type = character(), # Set at initialisation
    .numE = integer(), # Set at initialisation
    .d_time = numeric(), # Set at initialisation
    .transmission_type = character(), # Set at initialisation

    .time = 0,
    .beta = 0.05,
    .omega = 0.05,
    .gamma = 0.025,
    .delta = 0.005,
    .vacc = 0.001,
    .repl = 0.0001,
    .cull = 0.002,
    .trans_external = 0,
    .group_name = NA_character_,

    ## To store state:
    .saved_public = list(),
    .saved_private = list(),

    ## Private methods:
    check_state = function(){
      ## Use the new compartment_rule private utility method (see below):
      qassert(private$.S, private$compartment_rule())
      qassert(private$.E, private$compartment_rule(private$.numE))
      qassert(private$.I, private$compartment_rule())
      qassert(private$.R, private$compartment_rule())
      qassert(private$.N, private$compartment_rule())

      calcN <- private$.S + sum(private$.E) + private$.I + private$.R
      if(private$.update_type=="stochastic"){
        stopifnot(calcN == private$.N)
      }else if(private$.update_type=="deterministic"){
        stopifnot(all.equal(calcN, private$.N))
      }else{
        stop("Internal logic error")
      }
    },

    get_transmission_rate = function(){
      if(self$transmission_type=="frequency"){
        trans_internal <- private$.beta * private$.I / private$.N
      }else if(self$transmission_type=="density"){
        trans_internal <- private$.beta * private$.I
      }else{
        stop("Unrecognised transmission type: ", self$transmission_type)
      }
      trans_internal + self$trans_external
    },

    reset_N = function(){
      private$.N <- private$.S + sum(private$.E) + private$.I + private$.R
      ## Run the check_state method to ensure everything is OK:
      private$check_state()
    },

    ## A utility function to get the correct rule for checking S/E/I/R/N:
    compartment_rule = function(length=1){
      qassert(length, "X1(0,)")
      if(private$.update_type=="stochastic"){
        tt <- "X"
      }else if(private$.update_type=="deterministic"){
        tt <- "N"
      }else{
        stop("Internal logic error")
      }
      str_c(tt, length, "[0,)")
    }
  ),

  ## Active binding functions:
  active = mlist(

    #' @field S number of susceptible animals
    S = function(value){
      if(missing(value)) return(private$.S)
      qassert(value, private$compartment_rule())
      private$.S <- value
      private$reset_N()
    },
    E = function(value){
      if(missing(value)) return(sum(private$.E))
      qassert(value, private$compartment_rule())
      ## For E we distribute input randomly/equally over sub-compartments:
      if(private$.update_type=="stochastic"){
        private$.E <- rmultinom(1, value, rep(1/private$.numE,private$.numE))[,1]
      }else if(private$.update_type=="deterministic"){
        private$.E <- rep(value/private$.numE, private$.numE)
      }else{
        stop("Internal logic error")
      }
      private$reset_N()
    },
    I = function(value){
      if(missing(value)) return(private$.I)
      qassert(value, private$compartment_rule())
      private$.I <- value
      private$reset_N()
    },
    R = function(value){
      if(missing(value)) return(private$.R)
      qassert(value, private$compartment_rule())
      private$.R <- value
      private$reset_N()
    },

    ## Same as before
    beta = function(value){
      if(missing(value)) return(private$.beta)
      assert_number(value, lower=0)
      private$.beta <- value
    },
    omega = function(value){
      if(missing(value)) return(private$.omega)
      assert_number(value, lower=0)
      private$.omega <- value
    },
    gamma = function(value){
      if(missing(value)) return(private$.gamma)
      assert_number(value, lower=0)
      private$.gamma <- value
    },
    delta = function(value){
      if(missing(value)) return(private$.delta)
      assert_number(value, lower=0)
      private$.delta <- value
    },
    vacc = function(value){
      if(missing(value)) return(private$.vacc)
      assert_number(value, lower=0)
      private$.vacc <- value
    },
    repl = function(value){
      if(missing(value)) return(private$.repl)
      assert_number(value, lower=0)
      private$.repl <- value
    },
    cull = function(value){
      if(missing(value)) return(private$.cull)
      assert_number(value, lower=0)
      private$.cull <- value
    },
    time = function(){
      private$.time
    },

    ## N is now stored internally, but is still read only:
    N = function(){
      private$.N
    },

    ## As before, but including E:
    state = function(){
      tibble(Time = self$time, S = self$S, E = self$E, I = self$I, R = self$R)
    },

    ## New active binding functions:
    trans_external = function(value){
      if(missing(value)) return(private$.trans_external)
      qassert(value, "N1[0,)")
      private$.trans_external <- value
    },

    transmission_type = function(value){
      if(missing(value)) return(private$.transmission_type)
      private$.transmission_type <- match.arg(value,
        choices=c("frequency","density"))
    }
  )
)
